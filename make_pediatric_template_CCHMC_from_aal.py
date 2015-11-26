# -*- coding: utf-8 -*-
"""
Created on Tue Sep 29 14:51:35 2015

@author: edogerde
"""

# System import
import os
import subprocess


"" "make specific atlas base on the CCHMC of differents age and the AAL atlas"""

"""register AAL template on CCHMC specifics templates """
def template_to_child_template (template_adult,atlas_adult,template_children,
                                  path):
# Output directory
    trans_aff = os.path.join(path, "transf_aff.txt")
    template_adult_transaff = os.path.join(path, "template_adult_transaff.nii")
    template_adult_transnl = os.path.join(path, "template_adult_image_trans_nl.nii")
    trans_nl = os.path.join(path, "template_adult_trans_nl.nii")
    
    cmd = ["flirt", "-cost", "normmi", "-omat", trans_aff, "-in", template_children,
               "-ref",template_adult , "-out", template_adult_transaff]
    print "Executing: '{0}'.".format(" ".join(cmd))
    subprocess.check_call(cmd)
    
    # NL registration
    cmd = ["fnirt", "--ref={0}".format(template_adult),
           "--in={0}".format(template_children),
           "--iout={0}".format(template_adult_transnl),
           "--fout={0}".format(trans_nl), "--aff={0}".format(trans_aff),
           "--config=T1_2_MNI152_2mm.cnf"]
    # make my own configuration
    #"--config=/neurospin/grip/protocols/MRI/dosimetry_elodie_2015/atlas/atlas_aal/elodie.cnf"] 
    print "Executing: '{0}'.".format(" ".join(cmd))
    subprocess.check_call(cmd)
    
    return trans_aff, trans_nl

"""apply the matrix tranformation of the atlas on the CCHMC"""
#def atlas_to_child_template (atlas_adult, template_children, trans_nl, path):

# Output directory
    #registered_atlas = os.path.join(path, "registered_atlas.nii")
    #invt1_trans = os.path.join(path, "invt1_trans.nii")
    # Invert the nl warp
    #cmd = ["invwarp", "-w", trans_nl, "-o", invt1_trans, "-r",template_children ]
    #print "Executing: '{0}'.".format(" ".join(cmd))
    #subprocess.check_call(cmd)
    
    # Warp the labels
    #cmd = ["applywarp", "-i",atlas_adult , "-o", registered_atlas,
         #"-r",template_children , "-w",trans_nl , "--interp=nn"]
    #print "Executing: '{0}'.".format(" ".join(cmd))
   # subprocess.check_call(cmd)
    
    #return invt1_trans, registered_atlas
    #return registered_atlas

if __name__ == "__main__":
# Global parameters
    path = "/neurospin/grip/protocols/MRI/dosimetry_elodie_2015/atlas/atlas_aal"
    template_adult = os.path.join(path, "atlas_t1.nii.gz")
    atlas_adult = os.path.join(path, "atlas_labels.nii.gz")
    template_children = os.path.join(path,"CCHMC2_y.img")


trans_aff, trans_nl =  template_to_child_template(template_adult, atlas_adult,
                        template_children, path)

