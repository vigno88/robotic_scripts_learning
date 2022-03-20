import sys

# Tell the degree of freedom of possible joints
dof_joint = {'R':1,'P':1,'H':1,'C':2,'U':2,'S':3}

# Get the type of mechanism
# Returns true if it is spatial, false if planar
def get_mechanism_type(lines):
  for line in lines:
    line = line.strip()
    if line.startswith('.') and line[1:] == "planar":
      return False
    if line.startswith('.') and line[1:] == "spatial":
      return True
  print("Please specify the mechanism type at the top of the config file")
  quit()

# Parses the lines that are passed to it as an array
# Returns a list of links and a list of joints
def parse_file(lines):
  links = []
  joints = []  
  for i,line in enumerate(lines):
    line = line.strip()
    # Check if line is comment of mechanism type
    if line.startswith('#') or line.startswith('.') :
      continue
    line = line.replace(" ", "")
    # Check for empty line
    if len(line) == 0:
      continue

    tokens = line.split('-')
    if len(tokens) != 3: # Check if lines contains exactly 3 parts
      print("Line %d: all lines should be composed of two links joined by a joint" % i)
      quit()
    if tokens[1] not in ['R','P','H','C','U','S']: # Check if valid joint type
      print("Line %d: this joint %s is not valid" % (i,tokens[1]))
      quit()
    joints.append(tokens[1])
    if tokens[0] not in links: # Add first link
      links.append(tokens[0])
    if tokens[2] not in links: # Add second link
      links.append(tokens[2])
  return links,joints

# Takes a list of links, a list of joints and a flag that tell if the mechanism is
# spatial
# Returns the number of degrees of freedom
def compute_grubler(links, joints, is_spatial):
  sum_joints_freedom = 0
  for joint in joints:
    if not is_spatial and joint in ['H','C','U','S']:
      print("This type of joint: %s is not valid in planar config" % joint)
      quit()
    sum_joints_freedom += dof_joint[joint]
  if is_spatial:
    return 6*(len(links)-1-len(joints)) + sum_joints_freedom
  return 3*(len(links)-1-len(joints)) + sum_joints_freedom
      

def main():
  if len(sys.argv) != 2 or not sys.argv[1].endswith('.config'):
    print("The gruble calculator requires a single .config file as argument")
    quit()
  
  # Get the content of the file
  c_file = open(sys.argv[1],'r')
  lines = c_file.readlines()
  is_spatial = get_mechanism_type(lines)
  links, joints = parse_file(lines)
  
  # Compute the grubler formula
  dof = compute_grubler(links,joints,is_spatial)
  print("%s describes a mechanism with %d degrees of freedom" % (sys.argv[1],dof))

if __name__== "__main__":
  main()
