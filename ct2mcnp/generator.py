import os 
import numpy as np
import SimpleITK as sitk
from io import FileIO
from bisect import bisect_left
from scipy.ndimage import center_of_mass
from typing import List, Union


class CTVoxel:
    def __init__(self, ct:sitk.Image, material_config:dict[str, dict]) -> None:
        self.array = sitk.GetArrayFromImage(ct)
        self.spacing = np.around(np.asarray(ct.GetSpacing())/20, 3)
        self.size = np.asarray(ct.GetSize())
        
        self.hu_interval, self.element_table, self.density_table, self.keys = self._parse_material(material_config)
        self.array = np.clip(self.array, min(self.hu_interval), max(self.hu_interval))
        self.index_map = self._convert_hu_to_material()

    def _parse_material(self, material_config:dict[str, dict]):
        raw_hu_interval = list()
        element_table = dict()
        density_table = dict()
        keys = list()

        for i, material in material_config.items():
            raw_hu_interval += material['hu_interval']
            keys.append(int(i)) 
            element_table[int(i)] = dict(list(zip(material['nucleon'], material['fraction'])))
            density_table[int(i)] = material['density']
        hu_interval = list(set(raw_hu_interval))
        hu_interval.sort(key=raw_hu_interval.index)
        return hu_interval, element_table, density_table, keys
    
    def _convert_hu_to_material(self):
        index_map = np.zeros_like(self.array, dtype=np.uint32)
        it = np.nditer(self.array, flags=['multi_index'])
        while not it.finished:
            index = bisect_left(self.hu_interval, it[0])
            if index == 0: index = 1
            index_map[it.multi_index] = index
            it.iternext()
        return index_map


class MCNPGenerator:
    def __init__(self, ct:sitk.Image, config:dict, path) -> None:
        self.ct = ct 
        self.config = config 
        self.path = path 

    def run(self) -> None:
        self.init_MC_module()
        self.to_file()
        
    def init_MC_module(self):
        keys = self.config.keys()
        assert 'material' in keys and 'mode' in keys
        for card in ['source', 'tally', 'outcontrol']:
            if card not in keys:
                self.config[card] = None
        self.voxel = CTVoxel(self.ct, self.config['material'])
        self.material = Material(self.voxel)
        self.geom = Geometry(self.voxel, self.config['mode'])
        self.source = Source(self.config['source'])
        self.tally = Tally(self.geom, self.config['mode'], self.config['tally'])
        self.out = OutControl(self.config['outcontrol'])

    def to_file(self):
        with open(self.path, "w", encoding="utf-8") as f:
            f.write("c Geometry\n")
            self.geom.to_file(f)
            f.write("\n")
            f.write("c Data\n")
            self.material.to_file(f)
            self.source.to_file(f)
            self.tally.to_file(f)
            self.out.to_file(f)
            f.write("\n")


class Geometry:
    def __init__(self, voxel:CTVoxel, mode) -> None:
        self.mode = ','.join(mode)
        self.size = voxel.size
        self.spacing = voxel.spacing
        self.index_map = voxel.index_map
        self.keys = voxel.keys
        self.element_table = voxel.element_table
        self.density_table = voxel.density_table

        self._build_base_cell_surfaces()
        self._build_phantom_cell_surfaces()
        self._build_world_surface()

    def to_file(self, file:FileIO):
        for i in self.keys:
            file.write(f"{i}  {i}  {self.density_table[i]}")
            file.write(f"  {self.base_cell_geom}  u={i}  imp:{self.mode}=1\n")
            file.write(f"88{i}  0  #{i}  u={i}  imp:{self.mode}=1\n")
        file.write(f"998  0  {self.base_cell_geom} u=999 imp:{self.mode}=1\n")
        file.write(f"     lat=1  {self.fill}\n    ")
        accumulator = 1
        it = np.nditer(self.index_map, flags=['multi_index'])
        while not it.finished:
            file.write(f" {it[0]}")
            if accumulator % 20 == 0 and accumulator != np.prod(self.size):
                file.write("\n    ")
            elif accumulator % 20 == 0 and accumulator == np.prod(self.size):
                file.write("\n")
            elif accumulator == np.prod(self.size):
                file.write("\n")
            accumulator += 1
            it.iternext()
        file.write(f"999  0  {self.phantom_cell_geom}  fill=999 imp:{self.mode}=1\n")
        file.write(f"1000  1  -0.00129  -1000 #999 imp:{self.mode}=1\n")
        file.write(f"9999  0  1000  imp:{self.mode}=0\n")

        file.write('\n')
        file.writelines(self.base_cell_surfaces)
        file.writelines(self.phantom_cell_surfaces)
        file.write(self.word_cell_surface)

    def _build_base_cell_surfaces(self):
        self.base_cell_surfaces = (
            f"11    px    {self.spacing[0]:.3f}\n",
            f"12    px    {-self.spacing[0]:.3f}\n",
            f"13    py    {self.spacing[1]:.3f}\n",
            f"14    py    {-self.spacing[1]:.3f}\n",
            f"15    pz    {self.spacing[2]:.3f}\n",
            f"16    pz    {-self.spacing[2]:.3f}\n"
        )
        self.base_cell_geom = "-11 12 -13 14 -15 16"

    def _build_phantom_cell_surfaces(self):
        self.phantom_cell_surfaces = list()
        self.phantom_origin = list()
        self.phantom_limit = list()
        self.fill = "fill="
        half_size = self.size // 2
        if self.size[0] % 2 == 0:
            self.phantom_cell_surfaces.extend([
                f"111  px  {(self.size[0]+1)*self.spacing[0]:.3f}\n",
                f"112  px  {(-self.size[0]+1)*self.spacing[0]:.3f}\n"
            ])
            self.fill += f"{-half_size[0]+1}:{half_size[0]} "
            self.phantom_limit.append((self.size[0]+1)*self.spacing[0])
            self.phantom_origin.append((-self.size[0]+1)*self.spacing[0])
        else:
            self.phantom_cell_surfaces.extend([
                f"111  px  {self.size[0]*self.spacing[0]:.3f}\n",
                f"112  px  {-self.size[0]*self.spacing[0]:.3f}\n"
            ])
            self.fill += f"{half_size[0]}:{half_size[0]} "
            self.phantom_limit.append(self.size[0]*self.spacing[0])
            self.phantom_origin.append(-self.size[0]*self.spacing[0])
        if self.size[1] % 2 == 0:
            self.phantom_cell_surfaces.extend([
                f"113  py  {(self.size[1]+1)*self.spacing[1]:.3f}\n",
                f"114  py  {(-self.size[1]+1)*self.spacing[1]:.3f}\n"
            ])
            self.fill += f"{-half_size[1]+1}:{half_size[1]} "
            self.phantom_limit.append((self.size[1]+1)*self.spacing[1])
            self.phantom_origin.append((-self.size[1]+1)*self.spacing[1])
        else:
            self.phantom_cell_surfaces.extend([
                f"113  py  {self.size[1]*self.spacing[1]:.3f}\n",
                f"114  py  {-self.size[1]*self.spacing[1]:.3f}\n"
            ])
            self.fill += f"{half_size[1]}:{half_size[1]} "
            self.phantom_limit.append(self.size[1]*self.spacing[1])
            self.phantom_origin.append(-self.size[1]*self.spacing[1])
        if self.size[2] % 2 == 0:
            self.phantom_cell_surfaces.extend([
                f"115  pz  {(self.size[2]+1)*self.spacing[2]:.3f}\n",
                f"116  pz  {(-self.size[2]+1)*self.spacing[2]:.3f}\n"
            ])
            self.fill += f"{-half_size[2]+1}:{half_size[2]}"
            self.phantom_limit.append((self.size[2]+1)*self.spacing[2])
            self.phantom_origin.append((-self.size[2]+1)*self.spacing[2])
        else:
            self.phantom_cell_surfaces.extend([
                f"115  pz  {self.size[2]*self.spacing[2]:.3f}\n",
                f"116  pz  {-self.size[2]*self.spacing[2]:.3f}\n"
            ])
            self.fill += f"{half_size[2]}:{half_size[2]}"
            self.phantom_limit.append(self.size[2]*self.spacing[2])
            self.phantom_origin.append(-self.size[2]*self.spacing[2])
        self.phantom_cell_geom = "-111 112 -113 114 -115 116"

    def _build_world_surface(self):
        phantom_max_size = (self.size * self.spacing).max()
        world_radius = phantom_max_size + 100
        self.word_cell_surface = f"1000  so  {world_radius}\n"


class Material:
    def __init__(self, voxel:CTVoxel) -> None:
        self.element_table = voxel.element_table

    def to_file(self, file:FileIO):
        for index, elements in self.element_table.items():
            file.write(f"M{index}\n")
            file.write(self._parse_elements_dict(elements))
            
    def _parse_elements_dict(self, elements:dict):
        string = ''
        for nucleon, fraction in elements.items():
            string += f'     {nucleon} {fraction}\n'
        return string
    

class Source:
    def __init__(self, source_config:dict) -> None:
        self.source_str = 'SDEF\n'
        if source_config is None:
            return
        
        self.si_sp_str = ''
        self.si_sp_id = 0
        for key, value in source_config.items():
            if isinstance(value, dict):
                self.si_sp_id += 1
                self.source_str += f'     {key}=D{self.si_sp_id}\n'
                self.si_sp_str += self._parse_si_sp(value['si'], value['sp'])
            else:
                self.source_str += f'     {key}={_list_to_str(value)}\n'

    def to_file(self, file:FileIO):
        file.write(self.source_str)
        file.write(self.si_sp_str)

    def _parse_si_sp(self, si:list, sp:list):
        si_sp_str = f'#     SI{self.si_sp_id}     SP{self.si_sp_id}\n'
        for i, p in zip(si, sp):
            si_sp_str += f'     {i}     {p}\n'
        return si_sp_str


class Tally:
    def __init__(self, geom:Geometry, mode:list, tally_config:dict) -> None:
        self.config = tally_config
        if tally_config is None:
            return 
        
        self.mode_card = 'mode ' + ' '.join(mode)
        self.size = geom.size
        self.mesh_origin = geom.phantom_origin
        self.mesh_limit = geom.phantom_limit

    def to_file(self, file:FileIO):
        if self.config is None:
            return
        
        for id, value in self.config.items():
            file.write(f"fmesh{id}4:{value['particle']} geom=XYZ origin={self.mesh_origin[0]:.3f} {self.mesh_origin[1]:.3f} {self.mesh_origin[2]:.3f}\n")
            file.write(f"     imesh={self.mesh_limit[0]:.3f} iints={self.size[0]}\n")
            file.write(f"     jmesh={self.mesh_limit[1]:.3f} jints={self.size[1]}\n")
            file.write(f"     kmesh={self.mesh_limit[2]:.3f} kints={self.size[2]}\n")
            
            keys = value.keys()
            if 'de' in keys and 'df' in keys:
                file.write(self._parse_de_df(value['de'], value['df'], f'{id}4'))
            if 'fm' in keys:
                file.write(f'FM{id}4  {_list_to_str(value["fm"])}\n')
                
    def _parse_de_df(self, de, df, id):
        de_df_str = f'#     DE{id}     DF{id}\n'
        for e, f in zip(de, df):
            de_df_str += f'     {e}     {f}\n'
        return de_df_str


class OutControl:
    def __init__(self, control_config:dict) -> None:
        self.control_str = ''
        if control_config is None:
            return
        
        for key, value in control_config.items():
            self.control_str += f'{key} {_list_to_str(value)}\n'

    def to_file(self, file:FileIO):
        file.write(self.control_str) 
        
        
def _list_to_str(vec:list):
    if isinstance(vec, list):
        return ' '.join([str(i) for i in vec])
    else:
        return str(vec)
    
    
 
            