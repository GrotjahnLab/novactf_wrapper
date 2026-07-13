import os
import xml.etree.ElementTree as ET
import numpy as np
import glob

input_xml_dir = './'
tomo_dir = 'reconstruction/'
tomo_suffix = '.mrc_10.00Apx.mrc'
ref_name = "VLP_ref.mrc"
mask_name = "VLP_mask.mrc"
xmax = 1216
ymax = 1216
zmax = 318
sx = 2
sy = 2
sz = 1
gpu_list = "2,3"
angular_sampling = "11.00"
pixel_size = 10.00
output_dir = ref_name.split('.')[0] + '_' + angular_sampling + 'deg_output'
#make sure output directory is made
if not os.path.exists(output_dir):
    os.makedirs(output_dir)

xml_list = glob.glob(f'{input_xml_dir}/*xml')
with open("pytom.sbatch", 'w') as sbatch_file:
    for file_path in xml_list:
        print(file_path)
        #input file
        tree = ET.parse(file_path)
        root = tree.getroot()
        #getting relevant data
        dose = root.find('Dose').text.split('\n')
        tlt = root.find('Angles').text.split('\n')
        defocus = root.find('GridCTF')

        #converting text into arrays for defocus
        defocus_values = []
        for node in defocus.findall('Node'):
            defocus_values.append(np.int32(np.float64(node.get('Value'))*10000))
        defocus_values = np.array(defocus_values, dtype=float)
        print(defocus_values)

        #converting text into arrays for tlt and dose
        tlt_values = []
        for i in range(len(tlt)):
            if tlt[i] != '':
                tlt_values.append(np.float64(tlt[i]))
        tlt_values = np.array(tlt_values, dtype=float)
        dose_values = []
        for i in range(len(dose)):
            if dose[i] != '':
                dose_values.append(np.float64(dose[i]))
        dose_values = np.array(dose_values, dtype=float)

        #make the dose file name
        dose_file_name = file_path.replace('.xml', '.txt')
        #write the doses into a new txt file
        np.savetxt(dose_file_name, dose_values, fmt='%f')

        #make the tlt file name
        tlt_file_name = file_path.replace('.xml', '.tlt')
        #write the tlt into a new tlt file
        np.savetxt(tlt_file_name, tlt_values, fmt='%f')

        #write the defocus file
        defocus_file_name = file_path.split('.')[0] + '.defocus'
        with open(defocus_file_name, 'w') as f:
            for i in range(len(defocus_values)):
                if i ==0:
                    f.write(f"{i+1}\t{i+1}\t{tlt_values[i]}\t{tlt_values[i]}\t{defocus_values[i]}\t2")
                else:
                    f.write(f"\n{i+1}\t{i+1}\t{tlt_values[i]}\t{tlt_values[i]}\t{defocus_values[i]}")
        tomo_name = file_path.replace('.mrc.xml', tomo_suffix)
        sbatch_file.write(f"pytom_match_template.py  -v {tomo_dir}/{tomo_name}  -m {mask_name} -t {ref_name}  \
--angular-search {angular_sampling} --search-x 0 {xmax} --search-y 0 {ymax} --search-z 0 {zmax}  --voxel-size-angstrom  {pixel_size} \
--low-pass 25 --high-pass 1000 -d {output_dir } --per-tilt-weighting -a {tomo_dir}/{tlt_file_name}  \
--dose-accumulation {tomo_dir}/{dose_file_name} --defocus {tomo_dir}/{defocus_file_name} \
--voltage 300 --amplitude-contrast 0.1 --spherical-aberration 2.7 -g {gpu_list} -s {sx} {sy} {sz} \n")

    # Replace 'file.xml' with the path to your XML file
    #defocus_value, dose_value, tlt_value = parse_xml_file('JTL016A_3_L02_ts_003.mrc.xml')
    #create_output_files('JTL016A_3_L02_ts_003.mrc.xml', defocus_value, dose_value, tlt_value)
