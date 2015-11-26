import os
import subprocess

# Affine registration

atlas = "/neurospin/grip/protocols/MRI/dosimetry_elodie_2015/clemence/atlas/atlas_5_9/ANTS9-5Years3T_brain.nii.gz"
atlas_head = "/neurospin/grip/protocols/MRI/dosimetry_elodie_2015/clemence/atlas/atlas_5_9/ANTS9-5Years3T_head.nii.gz"
atlas_regis_head = "/neurospin/grip/protocols/MRI/dosimetry_elodie_2015/clemence/atlas/atlas_5_9/atlas_regis_head.nii"
trans_aff = "/neurospin/grip/protocols/MRI/dosimetry_elodie_2015/clemence/atlas/atlas_5_9/atlas_regis_head.txt"
label = "/neurospin/grip/protocols/MRI/dosimetry_elodie_2015/clemence/atlas/atlas_5_9/ANTS9-5Years3T_brain_ANTS_LPBA40_atlas.nii.gz"
label_regis_head = "/neurospin/grip/protocols/MRI/dosimetry_elodie_2015/clemence/atlas/atlas_5_9/label_regis_head.nii"
cmd = ["flirt", "-cost", "normmi", "-omat", trans_aff, "-in", atlas,
        "-ref", atlas_head, "-out", atlas_regis_head]
print "Executing: '{0}'.".format(" ".join(cmd))
subprocess.check_call(cmd)

cmd = ["flirt", "-in", label,
        "-ref", atlas_head, "-applyxfm", "-init", trans_aff, "-out", label_regis_head, "-interp", "nearestneighbour"]
print "Executing: '{0}'.".format(" ".join(cmd))
subprocess.check_call(cmd)
