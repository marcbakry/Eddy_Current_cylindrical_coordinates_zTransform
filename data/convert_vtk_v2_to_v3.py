# input file
f_in = open("./data/maillage.vtk")
# output file
f_out = open("./data/maillage_v3.vtk","w")

# header
f_in.readline()
f_out.write("# vtk DataFile Version 3.0\n")
# title of dataset
f_out.write(f_in.readline())
# ASCII or binary (force ASCII)
f_in.readline()
f_out.write("ASCII\n")
# Grid structure
f_in.readline()
f_out.write("DATASET UNSTRUCTURED_GRID\n")
# points
line  = f_in.readline()
sline = line.split(" ")
Npts  = int(sline[1])
f_out.write(line)
for i in range(0,Npts):
    f_out.write(f_in.readline())
# cells
line_cells = f_in.readline()
if not "CELLS" in line_cells:
    f_out.write(line_cells)
    line_cells = f_in.readline()
sline = line_cells.split(" ")
Ncells = int(sline[1])
cells = [line for line in [f_in.readline() for _ in range(Ncells)] if len(line)]
# cell types
line_cell_types = f_in.readline()
if not "CELL_TYPES" in line_cell_types:
    f_out.write(line_cell_types)
    line_cell_types = f_in.readline()
cell_types = [line for line in [f_in.readline() for _ in range(Ncells)] if len(line)]
# material ids
line_cell_ids = f_in.readline()
if not "CELL_DATA" in line_cell_ids:
    f_out.write(line_cell_ids)
    line_cell_ids = f_in.readline()
f_in.readline()
f_in.readline()
cell_ids = [line for line in [f_in.readline() for _ in range(Ncells)] if len(line)]
# count cells with reference greater than 2
type_min = 4
cells      = [cell for cell, cell_t in zip(cells,cell_types) if int(cell_t.strip())>=type_min]
cell_ids   = [cell_id for cell_id, cell_t in zip(cell_ids,cell_types) if int(cell_t.strip())>=type_min]
cell_types = [cell_t for cell_t in cell_types if int(cell_t.strip())>=type_min]
Ncells     = len(cells)
# total count
total = 0
for cell in cells:
    total += len(cell.split())

f_out.write("CELLS "+str(Ncells)+" "+str(total)+"\n")
for cell in cells:
    f_out.write(cell)

f_out.write("CELL_TYPES "+str(Ncells)+"\n")
for cell_t in cell_types:
    f_out.write(cell_t)

f_out.write("CELL_DATA "+str(Ncells)+"\n")
f_out.write("SCALARS MaterialID int 1\n")
f_out.write("LOOKUP_TABLE default\n")
for cell_id in cell_ids:
    f_out.write(cell_id)
f_out.write("\n")

# close
f_in.close()
f_out.close()