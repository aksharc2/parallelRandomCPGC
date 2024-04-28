import random
import sys

def save_to_mtx_file(nodes, p, file_name):
    rows = cols = nodes
    edges = 0
    print("\n Started generating graph")
    with open(file_name, 'w') as file:
        # Write header information
        file.write('%%MatrixMarket matrix coordinate real general\n')
        # Write matrix size information
        
        file.write(f'{rows} {cols} \n')
    
            
        for u in range(rows):
            for v in range(cols):
                if random.random() <= p/100: 
                    file.write(f'{u + 1} {v + 1} \n')
                    edges = edges + 1
    print("\n Graph generated")
    return edges


def edit_line_in_mtx(file_path, line_number, new_line):
    """
    Edit a specific line in a Matrix Market (.mtx) file.

    Parameters:
        file_path (str): The path to the .mtx file.
        line_number (int): The line number to be edited (1-based indexing).
        new_line (str): The new line content that will replace the old line.
    """
    print("\n Writing edges")
    # Read the original content from the file
    with open(file_path, 'r') as file:
        lines = file.readlines()

    # Check if the line number is within the range of the file
    if 1 <= line_number <= len(lines):
        # Update the specific line
        lines[line_number - 1] = new_line + '\n'

        # Rewrite the entire file with the updated content
        with open(file_path, 'w') as file:
            file.writelines(lines)
    else:
        print(f"Line number {line_number} is out of range.")
    print("\n Finished writing edges")

# Example usage:
if __name__ == '__main__':
    # Let's assume we have an existing .mtx file named "bipartite_matrix.mtx"
    if len(sys.argv) >= 2:
        nodes = int(sys.argv[1])
        p = int(sys.argv[2])
        exp= int(sys.argv[3])
        dataset_folder = 'datasets/'
        file_name = f'{dataset_folder}bipartite_graph_{nodes}_{p}_{exp}.mtx'
        edges = save_to_mtx_file(nodes, p, file_name)
        line_number = 2
        new_line = f'{nodes} {nodes} {edges}'
        edit_line_in_mtx(file_name, line_number, new_line)
    else:
        print("Please provide aruguments in the following order: nodes, density, experimentNo.")
