import os
import yaml

def create_crab_configs(samples_cfg, output_dir):
  os.makedirs(output_dir, exist_ok=True)
  with open(samples_cfg, 'r') as f:
    samples = yaml.safe_load(f)
  outputs = {}
  all_samples = set()
  era = None
  for sample_name, sample_desc in samples.items():
    try:
      if sample_name == 'GLOBAL':
        era = sample_desc['era']
        continue
      sample_type = sample_desc['sampleType']
      if sample_type not in outputs:
        outputs[sample_type] = {}
      outputs[sample_type][sample_name] = {
        'inputDataset': sample_desc['miniAOD'],
        'ignoreFiles': sample_desc.get('miniAOD_ignoreFiles', [])
      }
      all_samples.add(sample_desc['miniAOD'])
    except:
      print(f'Error while parsing {sample_name} entry.')
      raise
  for sample_type, sample_dict in outputs.items():
    mc_data = 'data' if sample_type == 'data' else 'mc'
    store_failed = sample_type != 'data'
    config = {
      'config': {
        'params': {
          'sampleType': mc_data,
          'era': era,
          'storeFailed': store_failed,
        }
      }
    }
    data = { }
    for name, desc in sample_dict.items():
      if len(desc['ignoreFiles']) == 0:
        data[name] = desc['inputDataset']
      else:
        data[name] = desc

    with open(os.path.join(output_dir, sample_type + '.yaml'), 'w') as f:
      yaml.safe_dump(config, f)
      f.write('\n')
      yaml.safe_dump(data, f)
  with open(os.path.join(output_dir, "all_samples.txt"), 'w') as f:
    for sample in sorted(all_samples):
      f.write(sample + '\n')

if __name__ == "__main__":
  import argparse
  parser = argparse.ArgumentParser(description='Create configs that will be used to submit crab jobs.')
  parser.add_argument('--samples', required=True, type=str, help="yaml file with sample list")
  parser.add_argument('--output', required=True, type=str, help="output path where crab configs will be stored")
  args = parser.parse_args()

  create_crab_configs(args.samples, args.output)

