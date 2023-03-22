import toml, argparse, os
import SimpleITK as sitk
from ct2mcnp.generator import CTVoxel, MCNPGenerator
   
if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument("-c", "--config", type=str, default="./config.toml", 
                        help="path to MCNP generator config file ended with 'toml'. Default: './config.toml'")
    parser.add_argument("-d", "--dirpath", type=str, default="./inp", help="path to the directory that contains output files. Default: './INP'")
    args = parser.parse_args()
    
    config_path = args.config
    output_path = args.dirpath
    with open(config_path, "r", encoding="utf-8") as f:
        config = toml.load(f)
        print("Have loaded MCNP generator config!")
        
    def generate_process(ct_path, config, out_path):
        ct_image = sitk.ReadImage(ct_path)
        generator = MCNPGenerator(ct_image, config, out_path)
        print(f"Generating {base_name} input file...", end="\t")
        generator.run()
        print("Complete")
    
    if not os.path.exists(output_path):
        os.makedirs(output_path)
    if isinstance(config["ct"], list):
        for one_ct_path in config['ct']:
            base_name = os.path.basename(one_ct_path).partition(".")[0]
            generate_process(one_ct_path, config, os.path.join(output_path, base_name))
            
    elif os.path.isdir(config['ct']):
        for ct_file in os.listdir(config['ct']):
            one_ct_path = os.path.join(config["ct"], ct_file)
            base_name = ct_file.partition(".")[0]
            generate_process(one_ct_path, config, os.path.join(output_path, base_name))
            
    elif isinstance(config["ct"], str):
        ct_path = config["ct"]
        base_name = os.path.basename(ct_path).partition(".")[0]
        generate_process(ct_path, config, os.path.join(output_path, base_name))
        
    else:
        raise ValueError("Invalid CT file paths!")