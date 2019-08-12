__author__ = 'phil'
import os
import struct
import numpy as np
from forestlas import lasIO_, lasStructure
import global_header

class recover_temp_files:

    def __init__(self, temp_directory, out_dir=False, ptFrmt=1, ):

        self.temp_directory = temp_directory
        if not out_dir:
            out_dir = self.temp_directory
        else:
            out_dir = out_dir
        
        if ptFrmt == 1:
            ptFrmt = lasStructure.point_fmt1struct()[0]
        dt = lasStructure.point_fmt1struct()[1]

        for las_temp in os.listdir(self.temp_directory):

            h = global_header.h()
            num_records = os.path.getsize(os.path.join(self.temp_directory, las_temp)) / h['pointreclen']
            rtn_num = [0, 0, 0, 0, 0]
            a = np.zeros(num_records, dtype=dt)

            with open(os.path.join(self.temp_directory, las_temp), 'rb') as fh:
                for j, line in enumerate(range(num_records)):
                    fh.seek(j * h['pointreclen'])
                    d = extract_return(ptFrmt, fh, h)
                    for field in d:
                        # print field
                        a[j][field] = d[field]

                    if 0 < d['rtn_num'] < 6:
                        rtn_num[d['rtn_num']-1] += 1

                    # break

            h['numptbyreturn'] = tuple(rtn_num)
            np2LAS(a, h, ptFrmt, outFile=os.path.join(out_dir, '{}.recovered.las'.format(las_temp[:-6])))

class recover_temp_file:

    def __init__(self, las_temp, out_dir=False, verbose=False):

        # self.temp_directory = temp_directory
        if not out_dir:
            out_dir, tile = os.path.split(las_temp)
        else:
            tile = os.path.split(las_temp)[1]

        ptFrmt = lasStructure.point_fmt1struct()[0]
        dt = lasStructure.point_fmt1struct()[1]
        h = global_header.h()
        num_records = os.path.getsize(las_temp) / h['pointreclen']
        rtn_num = [0, 0, 0, 0, 0]
        a = np.zeros(num_records, dtype=dt)

        with open(las_temp, 'rb') as fh:
            for j, line in enumerate(range(num_records)):
                fh.seek(j * h['pointreclen'])
                d = extract_return(ptFrmt, fh, h)
                for field in d:
                    # print field
                    a[j][field] = d[field]

                if 0 < d['rtn_num'] < 6:
                    rtn_num[d['rtn_num']-1] += 1

                # break

        h['numptbyreturn'] = tuple(rtn_num)
        np2LAS(a, h, ptFrmt,
               outFile=os.path.join(out_dir, '{}.recovered.las'.format(tile[:-5])),
               verbose=verbose)

def np2LAS(arr, h, ptFrmt, outFile=False, verbose=False):

    """
    Can be used to export a numpy array in .las format
    """

    h['gensoftware'] = 'CRC207 LiDAR analysis software  '
    h['sysid'] = 'CRC207 LiDAR analysis software  '
    h['xmin'] = arr['x'].min()
    h['xmax'] = arr['x'].max()
    h['ymin'] = arr['y'].min()
    h['ymax'] = arr['y'].max()
    h['zmin'] = arr['z'].min()
    h['zmax'] = arr['z'].max()
    ### sorting out the rtn_num tuple
    rtn = np.zeros(5).astype(int)
    for row in arr:
        rtn_num = row['rtn_num'] - 1
        if rtn_num < 5:
            rtn[rtn_num] += 1
    h['numptbyreturn'] = tuple(rtn)
    h['numptrecords'] = len(arr)
    h['guid2'] = 0

    with open(outFile, "wb") as out:
        for j in lasStructure.headerstruct():
            if j[2] == 'c':
                out.write(h[j[0]])
            elif j[3] > 1:
                out.write(struct.pack('=' + str(j[3]) + j[2], *h[j[0]]))
            else:
                out.write(struct.pack('=' + j[2] , h[j[0]]))

        ## write points
        out.seek(h['offset'])
        for d in arr:
            for i in ptFrmt:
                if i[0] == 'return_grp':
                    byte = ((d['scan_edge'] & 1) << 7) | ((d['scan_dir'] & 1) << 6) | ((d['rtn_tot'] & 7) << 3) | (d['rtn_num'] & 7)
                elif i[0] == 'x':
                    byte = (d['x'] - h['xoffset']) / h['xscale']
                elif i[0] == 'y':
                    byte = (d['y'] - h['yoffset']) / h['yscale']
                elif i[0] == 'z':
                    byte = (d['z'] - h['zoffset']) / h['zscale']
                else:
                    byte = d[i[0]]
                out.write(struct.pack('=' + i[2], byte))

    if verbose: print ".las file recoered to {}".format(outFile)

def extract_return(ptFrmt, fh, h):

    point_dictionary = {} # dictionary for storing values in for each point

    for ent in ptFrmt:
        byte = fh.read(ent[1])
        val = struct.unpack('=' + ent[2] , byte)[0]
        if ent[0] == 'x':
            val = (val * h['xscale'] ) + h['xoffset']
        if ent[0] == 'y':
            val = (val * h['yscale'] ) + h['yoffset']
        if ent[0] == 'z':
            val = (val * h['zscale'] ) + h['zoffset']
        if ent[0] == 'return_grp':
            point_dictionary['rtn_num'] = val & 7
            #if point_dictionary['rtn_num'] == 0:
            #    raise Exception
            point_dictionary['rtn_tot'] = (val >> 3) & 7
            point_dictionary['scan_dir'] = (val >> 6) & 1
            point_dictionary['scan_edge'] = (val >> 7)
            continue # required so that 'return_grp' is not added to dictionary

        point_dictionary[ent[0]] = val

    if point_dictionary["z"] > 1000000:
        raise NameError("z very high: {}".format(point_dictionary["z"]))

    return point_dictionary

if __name__ == '__main__':

    import multiprocessing
    import glob

    F = glob.glob('/var/folders/3g/x_8fg5zj7lvghkk1q3xz9d380000gp/T/lidar.processing.1402121.tmp/*.temp')
    # F = ((f, {'verbose':False}) for f in F)
    kw = {'verbose':True, 'out_dir':'/User/phil/ALS/SCRATCH'}

    # p = multiprocessing.Pool(4)
    # [p.apply(recover_temp_file, (f,), kw) for f in F]

    p = multiprocessing.Pool(4)
    p.map(recover_temp_file, ())
    for f in F:
        recover_temp_file(f[0], **f[1])
        # recover_temp_file(*arg, **kargs)