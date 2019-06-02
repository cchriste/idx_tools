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
fields=None
max_field_len=0

#****************************************************
def updateAndPublish(idx,field_idx,output,force=False):
    global fields
    field=fields[field_idx]
    if not force:
        global old_fields
        old_field=old_fields[field_idx]
        if (not old_field.minmax_set and field.minmax_set) or old_field != field:
            old_fields[field_idx]=copy.copy(field)
        else:
            return

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
    global I
    print(indent(I)+field_str)
        
#****************************************************
def setDefaultMinMax(field,dtype_str):
    # Can *almost* get normal limits from OpenVisus python wrappers,
    #   but not quite, so have to generate them myself.
    if field.ncomponents > 1:
        dtype_str=dtype_str[0:dtype_str.find('[')]

    from numpy import finfo,iinfo
    type_info = None
    if "int" in dtype_str:
        type_info = iinfo(dtype_str)
    else:
        type_info = finfo(dtype_str)

    if field.global_min is None:
        field.global_min = type_info.min * field.ncomponents
    if field.global_max is None:
        field.global_max = type_info.max * field.ncomponents

    print("field:",field.name,"dtype:",dtype_str,"range:","["+
          str(field.global_min)+",",str(field.global_max)+"]")


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
                field.min=field.max=val.item()
                field.minmax_set=True
                field.calc_minmax_func=calc_minmax
                # print("woo-hoo! min:",field.min,"max:",field.max)
                # print("wtf! global_min:",field.global_min,"global_max:",field.global_max)
            else:
                if val<field.min:
                    field.min=val.item()
                elif val>field.max:
                    field.max=val.item()

#****************************************************
# gets the largest box of the idx for this field that fits within the specified max memory size
# This version uses level box to properly compute bounds
def getLargestLevelBox(idx,field,lvl,maxsize):
    lbox=idx.getLevelBox(lvl)
    print("lvl:",lvl,"box:",lbox.toString(),end=', ')
    tooBig=True
    while tooBig:
        query=Query(idx,ord('r'))
        query.position=Position(lbox)
        query.start_resolution=lvl
        query.end_resolutions.clear()
        query.end_resolutions.append(lvl)
        idx.beginQuery(query)
        if query.getByteSize() > maxsize:
            lbox=NdBox(lbox.p1,lbox.middle())
        else:
            global max_field_len
            print("largest box for field %s at lvl %d is <%s> (%d bytes), overall region: <%s>, actual box: <%s>" %
                  (field.name.ljust(max_field_len),lvl,query.nsamples.toString(),query.getByteSize(),
                   query.position.getNdBox().size().toString(),query.position.getNdBox().toString()))
            tooBig=False
    #print("...and lbox itself (should be same):",lbox.toString())
    return lbox

#****************************************************
# gets the largest box of the idx for this field that fits within the specified max memory size
# *super* naive implementation. Should use pow2 dims or otherwise efficiently make idx queries
def getLargestBoxNaive(idx,field,lvl,maxsize):
    box=idx.getBox()
    tooBig=True
    while tooBig:
        query=Query(idx,ord('r'))
        query.position=Position(box)
        query.start_resolution=lvl
        query.end_resolutions.clear()
        query.end_resolutions.append(lvl)
        idx.beginQuery(query)
        if query.getByteSize() > maxsize:
            box=NdBox(box.p1,box.middle())
        else:
            global max_field_len
            print("largest box for field %s at lvl %d is <%s> (%d bytes), overall region: <%s>" %
                  (field.name.ljust(max_field_len),lvl,query.nsamples.toString(),query.getByteSize(),query.position.getNdBox().size().toString()))
            tooBig=False
    return box

#****************************************************
# gets the next box, returns None if there aren't any more
# TODO: [] try using a z-order to achieve faster convergence (measure and graph)
def getNextBox(idx,lvl,curr_box,box_size):
    idx_box=idx.getBox()
    next_box=NdBox(curr_box.p1,curr_box.p2)

    global I
    #print(indent(I)+"getNextBox(): curr_box:",curr_box.toString()+",","box_size:",box_size.toString())
    # x
    if next_box.p2[0] < idx_box.p2[0]:
        next_box.p1.set(0,     next_box.p1[0] + box_size[0])                 # increment x p1
        next_box.p2.set(0, min(next_box.p2[0] + box_size[0], idx_box.p2[0])) # increment x p2
        #print(indent(I)+"increment x, next_box:",next_box.toString())
        return next_box

    # y
    if next_box.p2[1] < idx_box.p2[1]:
        next_box.p1.set(1,     next_box.p1[1] + box_size[1])                 # increment y p1
        next_box.p2.set(1, min(next_box.p2[1] + box_size[1], idx_box.p2[1])) # increment y p2
        next_box.p1.set(0, idx_box.p1[0])                                    # wrap x p1
        next_box.p2.set(0, idx_box.p1[0] + box_size[0])                      # wrap x p2
        #print(indent(I)+"increment y, wrap x, next_box:",next_box.toString())
        return next_box
    
    # z
    if next_box.p2[2] < idx_box.p2[2]:
        next_box.p1.set(2,     next_box.p1[2] + box_size[2])                 # increment z p1
        next_box.p2.set(2, min(next_box.p2[2] + box_size[2], idx_box.p2[2])) # increment z p2
        next_box.p1.set(1, idx_box.p1[1])                                    # wrap y p1
        next_box.p2.set(1, idx_box.p1[1] + box_size[1])                      # wrap y p2
        next_box.p1.set(0, idx_box.p1[0])                                    # wrap x p1
        next_box.p2.set(0, idx_box.p1[0] + box_size[0])                      # wrap x p2
        #print(indent(I)+"increment z, wrap x and y, next_box:",next_box.toString())
        return next_box

    #print(indent(I)+"Done! No more boxes to read for this field at this level")
    return None


#****************************************************
I=0  # global indent level
def indent(I):
    return "%s"%('  '*I)

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
    global fields
    fields=[]
    global max_field_len
    #print("FIELDS TO CALCULATE:",fields_to_calculate)
    for field in idx.fields: 
        if not (fields_to_calculate is None or field.name in fields_to_calculate):
            #print("NOT calculating range for field",field.name)
            break

        field_item=SimpleNamespace(name=field.name,
                                   size=field.dtype.getByteSize(),
                                   dtype=field.dtype.toString(),
                                   ncomponents=field.dtype.ncomponents(),
                                   minmax_set=False,
                                   global_min=global_min,
                                   global_max=global_max,
                                   calc_minmax_func=calc_minmax_initial)
        setDefaultMinMax(field_item, field.dtype.toString())
        max_field_len=max(max_field_len,len(field.name))
        fields.append(field_item)
    # for field in fields:
    #     print("calculating range for:",field.name)

    global old_fields
    old_fields=copy.deepcopy(fields)
    
    # get box
    box=idx.getBox()

    # create access
    access=idx.createAccess()

    # get num resolution levels
    num_levels=idx.getMaxResolution()
    #num_levels=min(num_levels,8)   #debug: set a smaller max level (15 means 32^3 if dims >= 32^3)

    # get timesteps
    timesteps=idx.timesteps.getRange()
    #timesteps.To=timesteps.From+0  #debug: limit number of timesteps

    # tell 'em what we're gonna do
    print("\n** idxrange::calc_ranges **")
    print("\tdataset:",idxpath,"\n\tdims:",box.toString(),"\n\tnum_levels:",num_levels,"\n\ttimesteps:",timesteps.toString())
    print("\tusing maximum of",maxmem,"bytes during calculation")
    print("computing ranges for the following fields:")
    for field in fields:
        print("\tfield:",field.name,"dtype:",field.dtype)
        if field.ncomponents > 1:
            print("\t\t(ignoring because ncomponents is > 1)")  #todo: compute min/max for each component
            fields.remove(field)
            if len(fields) == 0:
                print("Nothing to compute!")
                return

    #debug: set getLargestBox functions
    getLargestBox=getLargestLevelBox    #this results in aligned queries
    #getLargestBox=getLargestBoxNaive   #this results in non-aligned queries and probably takes longer (todo: mesaure difference)
    idx_box=idx.getBox()
    
    # for each resolution level
    for lvl in range(0,num_levels+1):
        print("\nLevel:",lvl)

        # get box size for each field at this level
        largest_box_at_curr_level={}
        for field in fields:
            largest_box_at_curr_level[field.name]=getLargestBox(idx,field,lvl,maxmem)

        # for each timestep
        global I
        I += 1
        for ts in range(int(timesteps.From),
                        int(timesteps.To + timesteps.step),
                        int(timesteps.step)):
            print(indent(I)+"timestep:",ts)

            # read and compute range of all fields at this level
            I += 1 # increase indent
            for field in fields:
                if field.ncomponents > 1:
                    print("SKIPPING FIELD!!! (should decide to do this much earlier / todo)")
                    continue

                print(indent(I)+"field:",field.name)

                # get region size that fits in memory
                curr_box=largest_box_at_curr_level[field.name]
                box_size=curr_box.size()

                # make sure query box isn't bigger than idx box (otherwise we get default values!)
                curr_box.p2.set(2, min(curr_box.p2[2], idx_box.p2[2])) # increment z p2
                curr_box.p2.set(1, min(curr_box.p2[1], idx_box.p2[1])) # increment y p2
                curr_box.p2.set(0, min(curr_box.p2[0], idx_box.p2[0])) # increment x p2

                cnt=1
                I += 1 # increase indent
                while curr_box is not None:
#                    print(indent(I)+"Beginning query",str(cnt)+",","box:",curr_box.toString())
                    query=Query(idx,ord('r'))
                    query.field=idx.getFieldByName(field.name)
                    query.position=Position(curr_box)
                    query.start_resolution=lvl
                    query.end_resolutions.clear()
                    query.end_resolutions.append(lvl)
                    idx.beginQuery(query)
                    if query.canExecute():
                        print(indent(I)+"Executing query",str(cnt)+",","box:",curr_box.toString()+",","bounds:",query.nsamples.toString()+",","size:",query.getByteSize(),end="...")
                        idx.executeQuery(access,query)
                        #print("DONE executing query!")
                        if query.failed():
                            print(" failed! :-(");
                        else:
                            print(" succeeded! :-)");

                        #print("reading query buffer...")
                        buf=VisusKernelPy.Array.toNumPy(query.buffer)

                        # perform the comparisons
                        #print("performing minmax computation...")
                        field.calc_minmax_func(field,buf)
                    else:
                        # it's okay. Probably there just aren't any samples in this region at this level
                        print(indent(I)+"SKIPPING query",cnt,"(box:",curr_box.toString()+")","because",query.getLastErrorMsg())

                    # get next box (None if we're done)
                    #print("and getting the next box...")
                    curr_box=getNextBox(idx,lvl,curr_box,box_size)
                    #print("...which is",curr_box)

                    # publish current values
                    updateAndPublish(idx,fields.index(field),output)
                    cnt += 1
                    
                I -= 1 # reduce indent (queries)

            I -= 1 # reduce indent (field)
        I -= 1 # reduce indent (timestep)

    # publish all values
    print("\nFINISHED!\n")
    for i in range(0,len(fields)):
        updateAndPublish(idx,i,output,force=True)
    print("")

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
    parser.add_argument("--maxmem",required=False,type=int,default=2*1024*1024,help="maximum amount of memory to be used by this utility")
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
        

