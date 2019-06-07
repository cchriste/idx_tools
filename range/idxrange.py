#
# idxrange.py
#
# Incrementally computes range of fields within an idx volume.
#
#****************************************************

from OpenVisus import *
from numpy import nditer
import copy

class IdxRange:
    """Incrementally computes the range of a given local or remote IDX volume.
    Allows user to specify the specific field and timesteps for which to compute range, 
    as well as [min,max] values in which the computation will be limited."""

    def __init__(self, idx_path):
        # initialize instance variables
        #self.old_fields    = None
        self.changed       = True
        self.fields        = None
        self.max_field_len = 0
        self.path          = idx_path  # ISSUE: Dataset.url should just be exposed by SWIG
        self.I             = 0  # output indentation level (for debugging)

        # try to load IDX
        VisusIdxPy.IdxModule.attach()
        self.idx = LoadDataset(idx_path)
        if self.idx is None:
            raise ValueError("ERROR: cannot load %s" % idx_path)
        if not self.idx.valid():
            raise ValueError("ERROR: %s is invalid" % idx_path)

        self.idx_box=self.idx.getBox()

    def __del__(self):
        VisusIdxPy.IdxModule.detach()    

    def valid(self):
        return self.idx is not None
        
    def _indent(self):
        return "%s"%('  '*self.I)

    def updateAndPublish(self,field_idx,force=False):
        field=self.fields[field_idx]
        # if not force:
        #     old_field=self.old_fields[field_idx]
        #     if (not old_field.minmax_set and field.minmax_set) or old_field != field:
        #         self.old_fields[field_idx]=copy.copy(field)
        #     else:
        #         return

        # print the lines as they will appear in the .idx
        # example: SUANGSTR float32 compressed format(1) 
        field_str=field.name+" "+field.dtype
        try:
            field_str += " min("+str(field.min)+")"
        except:
            field_str += " min(<unset>)"
            None
        try:
            field_str += " max("+str(field.max)+")"
        except:
            field_str += " max(<unset>)"
        print(self._indent()+field_str)

    #****************************************************
    def setDefaultMinMax(self,field,dtype_str):
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
            field.global_min = type_info.min
        if field.global_max is None:
            field.global_max = type_info.max

        print("field:",field.name,"dtype:",dtype_str,"possible range:","["+
              str(field.global_min)+",",str(field.global_max)+"]")


    #****************************************************
    # def calc_minmax(self,field,data):
    def calc_minmax(self,field,data):
        for val in nditer(data):
            if val<field.min and val>=field.global_min:
                field.min=val.item()
            if val>field.max and val<=field.global_max:
                field.max=val.item()

    #****************************************************
    def calc_minmax_initial(self,field,data):
        for val in nditer(data):
            if val>=field.global_min and val<=field.global_max:
                if not field.minmax_set:
                    field.min=field.max=val.item()
                    field.minmax_set=True
                    field.calc_minmax_func=self.calc_minmax
                else:
                    if val<field.min:
                        field.min=val.item()
                    elif val>field.max:
                        field.max=val.item()

    #****************************************************
    # gets the largest box of the idx for this field that fits within the specified max memory size
    def getLargestLevelBoxDiv2(self,field,lvl,maxsize,curr_dim):
        lbox=self.idx.getLevelBox(lvl)
        bitmask=self.idx.getBitmask()
        idx_dim=self.idx.getPointDim()
        #print("lvl:",lvl,"box:",lbox.toString(),"idx_dim:",idx_dim)#,end=', ')

        if lvl > 0:
            curr_dim=bitmask[lvl]
        else:
            curr_dim=0              # example bitmask: V001010101012012012012012  (24 levels)

        tooBig=True
        i=1
        while tooBig:
            #print("curr_dim:",curr_dim,end=', ')
            query=Query(self.idx,ord('r'))
            query.position=Position(lbox)
            query.start_resolution=lvl
            query.end_resolutions.clear()
            query.end_resolutions.append(lvl)
            self.idx.beginQuery(query)
            #print("lbox:",lbox.toString(),"nsamples box:",query.nsamples.toString())
            if query.getByteSize() > maxsize:
                new_p2 = lbox.p2
                new_p2.set(curr_dim, int((lbox.p1[curr_dim] + lbox.p2[curr_dim])/2))
                lbox=NdBox(lbox.p1, new_p2)
                curr_dim = bitmask[lvl-i]
                i+=1
            else:
                # print("largest box for field %s at lvl %d is <%s> (%d bytes), overall region: <%s>, actual box: <%s>" %
                #       (field.name.ljust(self.max_field_len),lvl,query.nsamples.toString(),query.getByteSize(),
                #        query.position.getNdBox().size().toString(),query.position.getNdBox().toString()))
                tooBig=False
        return curr_dim,lbox

    #****************************************************
    # gets the next box, returns None if there aren't any more
    # TODO: [] try using a meta-z-order for the set of boxes of this size at this level to achieve even faster convergence (measure and graph)
    def getNextBox(self,lvl,curr_box,box_size):
        next_box=NdBox(curr_box.p1,curr_box.p2)

        #print(self._indent()+"getNextBox(): curr_box:",curr_box.toString()+",","box_size:",box_size.toString())
        # x
        if next_box.p2[0] < self.idx_box.p2[0]:
            next_box.p1.set(0,     next_box.p1[0] + box_size[0])                 # increment x p1
            next_box.p2.set(0, min(next_box.p2[0] + box_size[0], self.idx_box.p2[0])) # increment x p2
            #print(self._indent()+"increment x, next_box:",next_box.toString())
            return next_box

        # y
        if next_box.p2[1] < self.idx_box.p2[1]:
            next_box.p1.set(1,     next_box.p1[1] + box_size[1])                 # increment y p1
            next_box.p2.set(1, min(next_box.p2[1] + box_size[1], self.idx_box.p2[1])) # increment y p2
            next_box.p1.set(0, self.idx_box.p1[0])                                    # wrap x p1
            next_box.p2.set(0, min(self.idx_box.p1[0] + box_size[0], self.idx_box.p2[0]))  # wrap x p2
            #print(self._indent()+"increment y, wrap x, next_box:",next_box.toString())
            return next_box

        # z
        if next_box.p2[2] < self.idx_box.p2[2]:
            next_box.p1.set(2,     next_box.p1[2] + box_size[2])                 # increment z p1
            next_box.p2.set(2, min(next_box.p2[2] + box_size[2], self.idx_box.p2[2])) # increment z p2
            next_box.p1.set(1, self.idx_box.p1[1])                                    # wrap y p1
            next_box.p2.set(1, min(self.idx_box.p1[1] + box_size[1], self.idx_box.p2[1]))  # wrap y p2
            next_box.p1.set(0, self.idx_box.p1[0])                                    # wrap x p1
            next_box.p2.set(0, min(self.idx_box.p1[0] + box_size[0], self.idx_box.p2[0]))  # wrap x p2
            #print(self._indent()+"increment z, wrap x and y, next_box:",next_box.toString())
            return next_box

        #print(self._indent()+"Done! No more boxes to read for this field at this level")
        return None


    #****************************************************
    def calc_ranges(self,fields_to_calculate,global_min,global_max,maxmem,timestep):

        from types import SimpleNamespace

        # get fields
        self.fields=[]
        #print("FIELDS TO CALCULATE:",fields_to_calculate)
        for field in self.idx.fields: 
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
                                       calc_minmax_func=self.calc_minmax_initial)
            self.setDefaultMinMax(field_item, field.dtype.toString())
            self.max_field_len=max(self.max_field_len,len(field.name))
            self.fields.append(field_item)
        # for field in self.fields:
        #     print("calculating range for:",field.name)

        #self.old_fields=copy.deepcopy(self.fields)

        # maxmem in megabytes
        maxmem = maxmem << 20

        # get box
        box=self.idx.getBox()

        # create access
        access=self.idx.createAccess()

        # get num resolution levels
        num_levels=self.idx.getMaxResolution()

        # get timesteps
        timesteps=self.idx.timesteps.getRange()
        if timestep is not None:
            if timestep < timesteps.From or timestep > timesteps.To:
                print("ERROR: specified timestep ("+str(timestep)+") is out of range ("+timesteps.From+", "+timesteps.To+")")
            else:
                timesteps.From=timesteps.To=timestep

        # DEBUG
        #maxmem=2048                    #debug: reduce maxmem to a tiny value
        #num_levels=min(num_levels,16)   #debug: set a smaller max level (15 means 32^3 if dims >= 32^3)
        getLargestBox=self.getLargestLevelBoxDiv2#these are aligned and 1/2 size of previous (so should be just like levels themselves)

        # tell 'em what we're gonna do
        print("\n** idxrange::calc_ranges **")
        print("\tdataset:",self.path,"\n\tdims:",box.toString(),"\n\tnum_levels:",num_levels,"\n\ttimesteps:",timesteps.toString())
        print("\tusing maximum of",maxmem,"bytes during calculation")
        print("computing ranges for the following fields:")
        for field in self.fields:
            print("\tfield:",field.name,"dtype:",field.dtype)
            if field.ncomponents > 1:
                print("\t\t(ignoring because ncomponents is > 1)")  #todo: compute min/max for each component
                self.fields.remove(field)
                if len(self.fields) == 0:
                    print("Nothing to compute!")
                    return

        curr_dim=0

        # for each resolution level
        for lvl in range(0,num_levels+1):
            print("\nLevel:",lvl)

            # get box size for each field at this level
            largest_box_at_curr_level={}
            for field in self.fields:
                curr_dim,largest_box_at_curr_level[field.name]=getLargestBox(field,lvl,maxmem,curr_dim)

            # for each timestep
            self.I += 1
            for ts in range(int(timesteps.From),
                            int(timesteps.To + timesteps.step),
                            int(timesteps.step)):
                print(self._indent()+"timestep:",ts)

                # read and compute range of all fields at this level
                self.I += 1 # increase indent
                for field in self.fields:
                    print(self._indent()+"field:",field.name)

                    # get region size that fits in memory
                    curr_box=largest_box_at_curr_level[field.name]
                    box_size=curr_box.size()

                    # make sure query box isn't bigger than idx box (otherwise we get default values!)
                    curr_box.p2.set(2, min(curr_box.p2[2], self.idx_box.p2[2])) # increment z p2
                    curr_box.p2.set(1, min(curr_box.p2[1], self.idx_box.p2[1])) # increment y p2
                    curr_box.p2.set(0, min(curr_box.p2[0], self.idx_box.p2[0])) # increment x p2

                    cnt=1
                    self.I += 1 # increase indent
                    while curr_box is not None:
                        print(self._indent()+"Beginning query",str(cnt)+",","box:",curr_box.toString())
                        query=Query(self.idx,ord('r'))
                        query.field=self.idx.getFieldByName(field.name)
                        query.position=Position(curr_box)
                        query.start_resolution=lvl
                        query.end_resolutions.clear()
                        query.end_resolutions.append(lvl)
                        self.idx.beginQuery(query)
                        if query.canExecute():
                            self.idx.executeQuery(access,query)
                            field.calc_minmax_func(field,VisusKernelPy.Array.toNumPy(query.buffer))
                        #else:
                            #print("ERROR: query",cnt,"failed because",query.getLastErrorMsg(),"(box:",curr_box.toString()+")")

                        # get next box (None if we're done)
                        curr_box=self.getNextBox(lvl,curr_box,box_size)
                        cnt += 1

                    # publish current values
                    self.updateAndPublish(self.fields.index(field))
                    self.I -= 1 # reduce indent (queries)

                self.I -= 1 # reduce indent (field)
            self.I -= 1 # reduce indent (timestep)

        # publish all values
        print("\nFINISHED!\n")
        for i in range(0,len(self.fields)):
            self.updateAndPublish(i,force=True)
        print("")


   
#****************************************************
if __name__ == '__main__':

    import argparse
    parser = argparse.ArgumentParser(description="Incrementally computes ranges of fields within an IDX volume.")
    parser.add_argument('idxpath',help="path/url of IDX volume")
    parser.add_argument("-f","--fields",required=False,nargs='+',help="list of fields for which to compute min/max (default: all fields)")
    minmax_help=" values will be ignored). NOTE: if different min/max required for each fields, please run this utility separated for each field"
    parser.add_argument("--min",dest='minval',type=float,required=False,help="minimum value to include (lower"+minmax_help)

    parser.add_argument("--max",dest='maxval',type=float,required=False,help="maxiumu value to include (higher"+minmax_help)
    parser.add_argument("--maxmem",required=False,type=int,default=2,help="maximum amount of memory (in MB) to be used by this utility")
    parser.add_argument("--timestep",required=False,type=int,help="timestep for which to compute range")
    parser.add_argument("--profile",required=False,action='store_true',help="profile execution time of this utility")
    args = parser.parse_args()

    if args.minval is not None and args.maxval is not None:
        if args.maxval < args.minval:
            parser.error("ERROR: --min must be <= --max, but "+str(args.minval)+" is > "+str(args.maxval))

    if args.maxmem < 1:
        parser.error("maxmem "+str(args.minval)+" must be > 0")

    # create the range calculator
    calculator = IdxRange(args.idxpath)
        
    # run the range calculator, profiling if requested
    execstr='calculator.calc_ranges(fields_to_calculate=args.fields, global_min=args.minval, global_max=args.maxval, maxmem=args.maxmem, timestep=args.timestep)'
    if args.profile:
        import cProfile
        cProfile.run(execstr)
    else:
        exec(execstr)
        

