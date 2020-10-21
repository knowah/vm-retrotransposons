import os
import subprocess

def jf_path_check(jf_path, input_path=None):
  if not os.path.exists(jf_path):
    raise FileNotFoundError("jf_path ({}) does not exist.".format(jf_path))
  if not os.access(jf_path, os.X_OK):
    raise PermissionError("jf_path ({}) is not executable.".format(jf_path))
  if input_path and not os.path.exists(input_path):
    raise FileNotFoundError("input_path ({}) does not exist.".format(input_path))
  

def run_jellyfish_bc(fasta_path, jf_path, K, outfile, threads=1, other_opts=[]):
  # verify jellyfish installed and executable
  jf_path_check(jf_path, fasta_path)
  
  # calculate the size parameter for jellyfish bc
  # size = sum of (fasta sequence length minus K-1)
  bc_size = 0
  fasta_entries = 0
  with open(fasta_path, 'r') as inf:
  	for line in inf:
  		if line.startswith(">"):
  			fasta_entries += 1
  		else:
  			bc_size += len(line.strip())
  bc_size -= fasta_entries * (K-1)
  
  # run jellyfish bc
  jf_command = [jf_path, 'bc', '-m', str(K), '-s', str(bc_size), '-t', str(threads), '-o', outfile]
  jf_command += other_opts
  jf_command += [fasta_path]
  
  jf_proc = subprocess.run(jf_command, stdout=subprocess.DEVNULL, stderr=subprocess.PIPE)
  
  if jf_proc.returncode != 0:
    raise RuntimeError("jellyfish bc failed: {}".format(jf_proc.stderr.decode("utf-8")))
    
def run_jellyfish_count(fasta_path, jf_path, K, initial_hash_size, outfile, threads=1, bcfile=None, other_opts=[]):
  # verify jellyfish installed and executable
  jf_path_check(jf_path, fasta_path)
  
  jf_command = [jf_path, 'count', '-m', str(K), '-s', str(initial_hash_size), '-t', str(threads), '-o', outfile]
  if bcfile:
  	jf_command += ['--bc', bcfile]
  jf_command += other_opts
  jf_command += [fasta_path]
  
  jf_proc = subprocess.run(jf_command, stdout=subprocess.DEVNULL, stderr=subprocess.PIPE)
  
  if jf_proc.returncode != 0:
    raise RuntimeError("jellyfish count failed: {}".format(jf_proc.stderr.decode("utf-8")))

def run_jellyfish_dump(count_file, jf_path, outfile, other_opts=['-c', '-t']):
	# verify jellyfish installed and executable
	jf_path_check(jf_path, count_file)
	
	jf_command = [jf_path, 'dump', '-o', outfile]
	jf_command += other_opts
	jf_command += [count_file]
	
	jf_proc = subprocess.run(jf_command, stdout=subprocess.DEVNULL, stderr=subprocess.PIPE)

	if jf_proc.returncode != 0:
		raise RuntimeError("jellyfish count failed: {}".format(jf_proc.stderr.decode("utf-8")))
