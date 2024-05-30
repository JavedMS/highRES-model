"""
Change import capacities in BASE_gen.dd file
"""


def import_input_change(importx):

    with open(snakemake.input[0], "r", encoding="utf8") as file:
        list_of_lines = file.readlines()
    
    # Modifying all lines containing 'import.UP/.FX'
    for line_number, line in enumerate(list_of_lines):
        if '.Import.UP' in line or '.Import.FX' in line:
            parts = line.split()
            # Ensure that the line is correctly formatted with an expected identifier
            if len(parts) >= 2 and (parts[0].endswith('.Import.UP') or parts[0].endswith('.Import.FX')):
                # Convert the current value to float and multiply by importx
                value = float(parts[-1])
                new_value = value * importx
                # Reconstruct the line with the new value
                list_of_lines[line_number] = f'{parts[0]} {new_value}\n'         

    with open(snakemake.output.BASEgennewdd, "w", encoding="utf8") as file:
        file.writelines(list_of_lines)
        print(f"File written to {snakemake.output.BASEgennewdd}")


import_input_change(snakemake.params.ImportX)
