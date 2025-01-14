import sys
import os
import os.path
import xml.etree.ElementTree as ET
from base64 import b64decode
from array import array
import matplotlib.pyplot as plt

# Check if correct number of command-line arguments are provided
if len(sys.argv) < 4:
    print('Provide mzXML file, scan number and peptide sequence')
    sys.exit(1)
    
XML_file = sys.argv[1]
scan_number = sys.argv[2]
peptide_seq = sys.argv[3]
peptide_seq = peptide_seq.upper()

# Check if the file exists
if not os.path.isfile(XML_file):
    print("Error: File not found.")
    sys.exit(1)

# Check for non-integer scan number
try:
    int(scan_number)
except ValueError:
    print('Invalid input. Scan number must be an integer')
    sys.exit(1)

    
ns = '{http://sashimi.sourceforge.net/schema/}'

mzs = []
ints = []

# Parse through file to access the peaks
scan_num_found = False
for event, ele in ET.iterparse(XML_file):
    if ele.tag == ns + 'scan':
        if ele.attrib.get('num') == scan_number:
            scan_num_found = True
            peaks_element = ele.find(ns + 'peaks')
            if peaks_element is not None:
                peaks = array('f', b64decode(peaks_element.text))
                if sys.byteorder != 'big':
                    peaks.byteswap()
                mzs = peaks[::2]
                ints = peaks[1::2]
            # Check peaks for scan number   
            else:
                print("Error: Peaks not found for scan number", scan_number)
                sys.exit(1)
        ele.clear()

#Check if the scan number was found
if not scan_num_found:
    print("Error: Scan number", scan_number, "not found in the XML file.")
    sys.exit(1)



#Dictionary for molecular weights. 
mw = {'A': 71.04, 'C': 103.04, 'D': 115.04, 'E': 129.04, 'F': 147.07,
      'G': 57.02, 'H': 137.06, 'I': 113.08, 'K': 128.09, 'L': 113.08,
      'M': 131.04, 'N': 114.04, 'P': 97.05, 'Q': 128.06, 'R': 156.10,
      'S': 87.03, 'T': 101.05, 'V': 99.07, 'W': 186.08, 'Y': 163.06 }

# Compute b and y ions
def compute_ions(peptide):
    b_ions = []
    y_ions = []
    N_terminal = 0
    C_terminal = 0

    for i in range(len(peptide)):
        try:
            N_terminal += mw[peptide[i]]  
            b_ions.append(N_terminal + 1)  

            C_terminal += mw[peptide[-(i + 1)]]  
            y_ions.append(C_terminal + 19)
        except KeyError:
            print('Peptide seqeunce contains an invalid amino acid.')
            sys.exit(1)
    return b_ions, y_ions

b_ions, y_ions = compute_ions(peptide_seq)

#print(b_ions,y_ions)


#Annotate b and y ions to match m/z values to peaks
def annotate_peaks(b_ions, y_ions, mzs, ints):
    b_annots = []
    y_annots = []

    intensity_threshold = max(ints) * 0.05
    tolerance = 0.05

    for mz, intensity in zip(mzs, ints):
        valid_peak = intensity >= intensity_threshold

        if valid_peak:
            for j, b_ion in enumerate(b_ions):
                if abs(mz - b_ion) <= 0.05: #Tolerance
                    annotation = ('b' + str(j + 1), mz, intensity)
                    b_annots.append(annotation)#matching b-ions

            for j, y_ion in enumerate(y_ions):
                if abs(mz - y_ion) <= 0.05:
                    annotation = ('y' + str(j + 1), mz, intensity)
                    y_annots.append(annotation)#matching y-ions

    return b_annots, y_annots

b_annots, y_annots = annotate_peaks(b_ions, y_ions, mzs, ints)

if not b_annots and not y_annots:
    print('No matching b or y ions found in the spectrum.', file=sys.stderr)
    sys.exit(1)
elif not b_annots:
    print('No matching b ions found in the spectrum.', file=sys.stderr)
    sys.exit(1)
elif not y_annots:
    print('No matching y ions found in the spectrum.', file=sys.stderr)
    sys.exit(1)

#print("b-ion annotations:", b_annots)#matched b-ions
#print("y-ion annotations:", y_annots)#matched y-ions



# Extract labels, m/z values, and intensities for b and y ions for plotting
b_labels = []
b_mzs = []
b_ints = []

for matched_b in b_annots:
    b_labels.append(matched_b[0])
    b_mzs.append(matched_b[1])
    b_ints.append(matched_b[2])

y_labels = []
y_mzs = []
y_ints = []

for matched_y in y_annots:
    y_labels.append(matched_y[0])
    y_mzs.append(matched_y[1])
    y_ints.append(matched_y[2])


#Plotting unmatched and matched peaks
plt.stem(mzs, ints, linefmt="black", markerfmt="None", label="Unmatched Peaks")
plt.stem(b_mzs, b_ints, linefmt="red", markerfmt="None", label="Matched b-ions")
plt.stem(y_mzs, y_ints, linefmt="blue", markerfmt="None", label="Matched y-ions")


# Plot annotations for b ions
for i in range(len(b_labels)):
    plt.text(b_mzs[i], b_ints[i], b_labels[i], ha="center", fontsize=8, color="red")

# Plot annotations for y ions
for i in range(len(y_labels)):
    plt.text(y_mzs[i], y_ints[i], y_labels[i], ha="center", fontsize=8, color="blue")

# Customize plot
plt.xlabel('m/z')
plt.ylabel('Intensity')
plt.title('Spectrum peaks for '+ 'Scan ' + scan_number + ' & Peptide: ' + peptide_seq)
plt.legend(['Unmatched Peaks', 'Matched b-ions', 'Matched y-ions'])

plt.show()









        


    




























