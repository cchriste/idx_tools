#
# idxrange.py
#
# Incrementally computes range of fields within an idx volume.
#
#****************************************************

from OpenVisus import *
from numpy import nditer
import copy

#****************************************************
old_fields=None

#****************************************************
def updateAndPublish(idx,fields,output):
    global old_fields
    same=True
    for old_field, field in zip(old_fields, fields):
        if (not old_field.minmax_set and field.minmax_set) or old_field != field:
            same=False
            break

    if same:
        return
    old_fields=copy.deepcopy(fields)

    # todo: maybe field.setParam('min',field.min)... not sure since there are no
    #       params at all in the idx field in the idx returned by LoadDataset.
    #
    # for field in fields: idx.field.minval=field.minval
    #     idx.field.maxval=field.maxval
    #     field.publish()

    #if output is a remote idx, server.updateDataset(output,idx)   (todo)
    #else:
    # idx.save()   #need a lock? (maybe somewhere)
    # Use something like this to just look at the idx file as it updates:
    # https://unix.stackexchange.com/questions/101271/open-a-text-file-in-a-terminal-and-auto-refresh-it-whenever-it-is-changed

    # for now, just print the values to stdout
    # todo: next step is to write them to separate files for each level, which will enable a graph showing rate of convergence
    #       (use this as evidence of non-linear progression)
    #
    # print the lines as they will appear in the .idx
    # example: SUANGSTR float32 compressed format(1) 
    for field in fields:
        if field.calc_range:
#            if fields.index(field) > 0:
#                print('+',end=' ')
            #todo: not sure how to get things like 'compressed' and 'format(0)' from Field, but gonna ignore it for now
            field_str=field.name+" "+field.dtype
            try:
                field_str += " min("+str(field.min)+")"
            except:
                # print("still no min in field:",field)
                # field_str += " min("+field.min+")"
                None
            try:
                field_str += " max("+str(field.max)+")"
            except:
                None
            print(field_str)
        
#****************************************************
# Returns tuple of the next box to be read at the requested level and the number read so far.
# There might need to be multiple reads to do the whole dataset in case it's bigger than maxsize.
# Use a meta z-order so the regions incrementally approach the final solution more quickly (todo).
# Returns 0 for nRead when done
#
def getBox(idx,field,lvl,maxsize,nRead): 
    #todo: compute size... there must be something that does this since it's just like a region query at a given resolution for the viewer
    box=dataset.getBox()
    nRead=0  # 0 means this is the last box and we're done reading this field at this level
    return box,nRead

#****************************************************
def setDefaultMinMax(field,dtype_str):
    # Can *almost* get normal limits from python wrappers, but not quite, so have to generate them myself.
    from numpy import finfo,iinfo
    type_info = None
    print("getting min/max for dtype",dtype_str)
    if field.ncomponents > 1:
        dtype_str=dtype_str[0:dtype_str.find('[')]
    if "int" in dtype_str:
        type_info = iinfo(dtype_str)
    else:
        type_info = finfo(dtype_str)
    if field.global_min is None:
        field.global_min = type_info.min
    if field.global_max is None:
        field.global_max = type_info.max

    print("global_min for",dtype_str,"is",field.global_min)
    print("global_max for",dtype_str,"is",field.global_max)


#****************************************************
def calc_minmax(field,data):
    for val in nditer(data):
        if val>=field.global_min and val<=field.global_max:
            if val<field.min:
                field.min=val.item()
            elif val>field.max:
                field.max=val.item()

#****************************************************
def calc_minmax_initial(field,data):
    for val in nditer(data):
        if val>=field.global_min and val<=field.global_max:
            if not field.minmax_set:
                #print("type(val):",type(val),val)
                field.min=field.max=val.item()
                field.minmax_set=True
                global_calc_minmax_func=calc_minmax
                #print("woo-hoo! min:",field.min,"max:",field.max)
            else:
                if val<field.min:
                    field.min=val.item()
                elif val>field.max:
                    field.max=val.item()

#****************************************************
global_calc_minmax_func=calc_minmax_initial

#****************************************************
def calc_ranges(idxpath,fields_to_calculate,global_min,global_max,maxmem,output):
    VisusIdxPy.IdxModule.attach()

    # read IDX
    idx=LoadDataset(idxpath)
    if idx is None:
        print("ERROR: cannot load %s" % idxpath)
        return -1
    if not idx.valid():
        print("ERROR: %s is invalid" % idxpath)
        return -1

    # initialize output volume
    # if len(output) > 0:
    #     try:
    #         outputidx=idx.open(output)
    #     except <...>:
    #         if isUrl(output):
    #             outputidx=createRemoteIdx(output)
    #         else:
    #             outputidx=idx.create(output)
    #     except:
    #         print("ERROR: cannot open %s" % output)
    # elif input is a writable file: (todo)
    #     make a backup copy of input idx
    #     outputidx=idx
    # elif input is a url and we have permission to modify: (todo)
    #     outputidx=idx

    from types import SimpleNamespace

    # get fields
    fields=[]
    for field in idx.fields: 
        field_item=SimpleNamespace(name=field.name,
                                   size=field.dtype.getByteSize(),
                                   dtype=field.dtype.toString(),
                                   ncomponents=field.dtype.ncomponents(),
                                   minmax_set=False,
                                   calc_range=fields_to_calculate is None,
                                   global_min=global_min,
                                   global_max=global_max)
        setDefaultMinMax(field_item, field.dtype.toString())
        if fields_to_calculate is not None:
            if field.name in fields_to_calculate:
                field_item.calc_range = True
        fields.append(field_item)
    
    global old_fields
    old_fields=copy.deepcopy(fields)
    
    # get box
    box=idx.getBox()

    # create access
    access=idx.createAccess()

    # get num resolution levels
    num_levels=idx.getMaxResolution()
    #debug: set a smaller max level (15 means 32^3 if dims >= 32^3)
    num_levels=min(num_levels,8)

    # get timesteps
    timesteps=idx.timesteps.getRange()

    # tell 'em what we're gonna do
    print(idxpath,"\n\tdims:",box.toString(),"\n\tnum_levels:",num_levels,"\n\ttimesteps:",timesteps.toString())
    print("computing ranges for the following fields:")
    for field in fields:
        if fields_to_calculate is not None and not field.name in fields_to_calculate:
            continue
        print("\tfield:",field.name,"dtype:",field.dtype)
        if field.ncomponents > 1:
            print("\t\t(ignoring because ncomponents is > 1)")

    # for each resolution level
    for lvl in range(0,num_levels+1):
        print("level:",lvl)

        # for each timestep
        for ts in range(int(timesteps.From),
                        int(timesteps.To + timesteps.step),
                        int(timesteps.step)):
            #print("timestep:",ts)

            # read and compute range of all fields at this level
            for field in fields:
                if field.ncomponents > 1:
                    continue
                #print("field:",field.name)

                query=Query(idx,ord('r'))
                query.position=Position(box)
                query.field=idx.getFieldByName(field.name)
                query.start_resolution=lvl
                query.end_resolutions.clear()
                query.end_resolutions.append(lvl)
                idx.beginQuery(query)
                # get region to read next (todo)
        #        nsamples=query.nsamples.innerProduct()
        #        size=field.size*nsamples
        #        if (too big...) todo
                idx.executeQuery(access,query)
                buf=VisusKernelPy.Array.toNumPy(query.buffer)

                # perform the comparisons
                global_calc_minmax_func(field,buf)
                # print(field)
                # print(field.min)
                # print(type(field.min))

            # publish current values
            updateAndPublish(idx,fields,output)

    # for each resolution level
    # for lvl in (1,num_levels):
    #     # read and compute range of all fields at this level
    #     for field in fields:
    #         while True:
    #             # get region to read next
    #             box,count=getBox(idx,field,lvl,maxsize,count=0)
    #             data=idx.read(field,box,fromh=lvl-1,toh=lvl)
    #             for i in data.size:
    #                 if i<field.minval and i>minval:
    #                     field.minval=i
    #                 else if i>field.maxval and i<maxval:
    #                     field.maxval=i
    #             # publish current values
    #             updateAndPublish(fields.ranges)

    #             # getBox returns 0 when it's the last box
    #             if not count:
    #                 break;


    VisusIdxPy.IdxModule.detach()    
 

   
#****************************************************
if __name__ == '__main__':

    import argparse
    parser = argparse.ArgumentParser(description="Incrementally computes ranges of fields within an IDX volume.")
    parser.add_argument('idxpath',help="path/url of IDX volume")
    parser.add_argument("-f","--fields",required=False,nargs='+',help="list of fields for which to compute min/max (default: all fields)")
    minmax_help=" values will be ignored). NOTE: if different min/max required for each fields, please run this utility separated for each field"
    parser.add_argument("--min",dest='minval',type=float,required=False,help="minimum value to include (lower"+minmax_help)

    parser.add_argument("--max",dest='maxval',type=float,required=False,help="maxiumu value to include (higher"+minmax_help)
    parser.add_argument("-o","--output",required=False,default="",help="path to output IDX file (default: same as input, or 'idx_with_ranges.idx' if input is not writable)")
    parser.add_argument("--maxmem",required=False,default=2*1024*1024,help="maximum amount of memory to be used by this utility")
    parser.add_argument("--profile",required=False,action='store_true',help="profile execution time of this utility")
    args = parser.parse_args()

    if args.minval is not None and args.maxval is not None:
        if args.maxval < args.minval:
            print("ERROR: --min must be <= --max, but",args.minval,"is >",args.maxval)
            exit(1)

    # run the range calculator, profiling if requested
    execstr='calc_ranges(idxpath=args.idxpath, fields_to_calculate=args.fields, global_min=args.minval, global_max=args.maxval, maxmem=args.maxmem, output=args.output)'
    if args.profile:
        import cProfile
        cProfile.run(execstr)
    else:
        exec(execstr)
        

