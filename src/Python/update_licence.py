#!/usr/bin/python

import os

def update_file(file_name, header, comment):
  #print file_name
  f = open(file_name, "r")
  old_lines = f.readlines()
  f.close()
  new_lines = []
  in_header = True
  for line in old_lines:
    if in_header:
      if not line.startswith(comment + " +"):
        in_header = False
    if not in_header:
      new_lines.append(line)
  f = open(file_name, "w")
  for line in header:
    f.write(comment + " " + line)
  for line in new_lines:
    f.write(line)
  

def recursive_traversal(diri, header):
  fns = os.listdir(dir)
  for fn in fns:
    file_name = os.path.join(dir,fn)
    if os.path.isdir(file_name):
      recursive_traversal(file_name)
    else:
      if (file_name.endswith(".h")):
        update_file(file_name, header, "//")
      if (file_name.endswith(".cpp")):
        update_file(file_name, header, "//")
      if (file_name.endswith(".cu")):
        update_file(file_name, header, "//")
      if (file_name.endswith(".pro")):
        update_file(file_name, header, "#")
      if (file_name.endswith(".sh")):
        update_file(file_name, header, "#")
      if (file_name.endswith(".bash")):
        update_file(file_name, header, "#")
      if (file_name.endswith(".py")):
        update_file(file_name, header, "#")


f = open("licence_header.txt")
header = f.readlines()
f.close

recursive_traversal('.')


