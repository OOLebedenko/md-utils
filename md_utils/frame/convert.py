import re

# rosetta hydrogen pattern to substitute
hydrogen_pattern = re.compile(r"\d[A-Z]+\d*")

def clean_rosetta_pdb(path_to_pdb):
    """
    cleans rosetta notation in pdb file to cinventional
    :param path_to_pdb: file to clean
    :return: writes .pdb with fixed atom notations and .txt with extra information from initial file
    """
    # extract structure name
    name = path_to_pdb.split('/')[-1][:-4]
    # define flags indicating changes
    separeted_data, fixed_atoms = False, False

    # read initial file
    with open(path_to_pdb, 'r') as f:
      lines = f.readlines()

    # separate structure and energy data
    pdb, data = [], []
    for i in range(len(lines)):
        if '# All' in lines[i]:
          sep_line = i
          separeted_data = True
          
    if separeted_data:
      pdb, data = lines[:sep_line], lines[sep_line:]
      with open(f'{name}_energy_data.txt', 'w') as f:
        f.writelines(data)
    else:
      pdb = lines
      print('Warning: no extra data was found in pdb')

    # clean hydrogen atom notation
    cleaned_pdb = []
    for line in pdb:
      result = hydrogen_pattern.search(line)
      if result is not None:
        h = result.group(0)[1:] + result.group(0)[0] # replace index
        line = re.sub(hydrogen_pattern, h, line) # substitute atom name
        fixed_atoms = True
      cleaned_pdb.append(line)

    if fixed_atoms:
      with open(f'{name}_cleaned.pdb', 'w') as f:
        f.writelines(cleaned_pdb)
    else:
      print('No hydrogen atoms requiring fixes were found')

if __name__ == "__main__":
    path_to_pdb = 'your path'
    clean_rosetta_pdb(path_to_pdb)
