# -*- coding: utf-8 -*-

"""
Resample the nii image (RD) from a voxel size xi,yi,zi to a voxel size xf, yf,zf
"""

import os
import glob


def convert_nii_to_minf(rd_path):
    sujet_minf = os.path.join(outputdir_Subj, "sujet_rd.nii.gz")
    file_minf = sujet_minf + '.minf'
    # AimsFileConvert
    cmd = "AimsFileConvert " + \
          " -i " + rd_path + \
          " -o " + sujet_minf
    
    os.system(cmd)
    return file_minf

 
def findsize(di, ri, rf):
    ''''''
    dfx = (ri[0]*di[0])/rf[0]
    dfy = (ri[1]*di[1])/rf[1]
    dfz = (ri[2]*di[2])/rf[2]
    df = [dfx, dfy, dfz]
    return df


def readminf(file_minf):
    globalVariables = dict()
    localVariables = dict()
    execfile(file_minf, globalVariables, localVariables)

    file_info = localVariables['attributes']
    resolution = file_info['voxel_size'][:-1]
    dimension = file_info['volume_dimension'][:-1]
    return resolution, dimension

def Resample(dim, res, rd_path, outputdir):
    # Output_file
    rd_resample = os.path.join(outputdir, "rd_resample.nii")
    # Resample
    cmd = "AimsResample " + \
          " -i " + rd_path + \
          " -o " + rd_resample + \
          " --dx " + str(dim[0]) + \
          " --dy " + str(dim[1]) + \
          " --dz " + str(dim[2]) + \
          " --sx " + str(res[0]) + \
          " --sy " + str(res[1]) + \
          " --sz " + str(res[2]) + \
          " -t " + str(0)
    print "Executing: %s" % (cmd)
    os.system(cmd)

    return rd_resample


if __name__ == "__main__":
    # Input File
    subject_path = "/neurospin/grip/protocols/MRI/dosimetry_elodie_2015/subject_24_rd_resample"
    rd_paths = glob.glob(os.path.join(subject_path, "sujet*", "rd", "*.nii.gz"))
    outputdir = os.path.join(subject_path, "results_rd+_resample_rd")
    if not os.path.exists(outputdir):
        os.makedirs(outputdir)
    resolution_final = [1, 1, 1]
            
    for rd_path in rd_paths:
        print "Processing {0}".format(rd_path)
        outputdir_Subj = os.path.join(outputdir, rd_path.split('/')[-3])
        if not os.path.exists(outputdir_Subj):
            os.makedirs(outputdir_Subj)
       
        sujet_minf = convert_nii_to_minf(rd_path)
        resolution, dimension = readminf(sujet_minf)
        dimension_finale = findsize(dimension, resolution, resolution_final)
        rd_resample = Resample(dimension_finale, resolution_final,
                               rd_path, outputdir_Subj)
