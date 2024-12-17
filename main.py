import pyvista as pv
from rdkit import Chem
from rdkit.Chem import AllChem
import os

# Crear carpeta de salida
output_folder = "imagenes_moleculas_pyvista"
os.makedirs(output_folder, exist_ok=True)

# Generar coordenadas 3D para RDKit
def get_3d_coords(mol):
    AllChem.EmbedMolecule(mol, AllChem.ETKDG())
    conf = mol.GetConformer()
    coords = [conf.GetAtomPosition(i) for i in range(mol.GetNumAtoms())]
    return [(atom.x, atom.y, atom.z) for atom in coords]

# Generar y guardar imágenes con PyVista
def generate_3d_images_with_pyvista(coords, output_folder, mol_name):
    for idx, angle in enumerate([(0, 0, 0), (90, 0, 0), (0, 90, 0), (45, 45, 0)]):
        plotter = pv.Plotter(off_screen=True)
        # Crear puntos y visualización
        cloud = pv.PolyData(coords)
        plotter.add_mesh(cloud, color="blue", point_size=20, render_points_as_spheres=True)
        plotter.camera_position = "iso"  # Vista isométrica inicial
        plotter.view_vector(angle)  # Rotar la vista

        # Guardar la imagen
        output_path = os.path.join(output_folder, f"{mol_name}_view_{idx+1}.png")
        plotter.screenshot(output_path)
        print(f"Imagen guardada: {output_path}")
        plotter.close()

# Procesar moléculas
sdf_file = "TAGORIZINE.sdf"
mol_supplier = Chem.SDMolSupplier(sdf_file)

for idx, mol in enumerate(mol_supplier):
    if mol:
        mol_name = f"mol_{idx+1}"
        mol_with_Hs = Chem.AddHs(mol)
        mol_output_folder = os.path.join(output_folder, mol_name)
        os.makedirs(mol_output_folder, exist_ok=True)
        
        coords = get_3d_coords(mol_with_Hs)
        generate_3d_images_with_pyvista(coords, mol_output_folder, mol_name)
