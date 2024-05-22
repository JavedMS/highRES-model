"""
Change transmission node capacities in the GAMS code
"""


def build_transmission(trans):
    """
    Crudely modify the GAMS code textfile, will break if line numbers change
    """
    with open(snakemake.input[0], "r", encoding="utf8") as file:
        list_of_lines = file.readlines()
    list_of_lines[72] = '$setglobal fx_trans "NO"'+"\n"

    if trans == "Default":
        with open(snakemake.output[0], "w", encoding="utf8") as file:
            file.writelines(list_of_lines)

    if trans == "FIXED":
        list_of_lines[328] = "\n" + 'var_trans_pcap.FX(z,z_alias,"HVAC400KV")$(trans_links(z,z_alias,trans)) = trans_links_cap(z,z_alias,trans);'+"\n"
        with open(snakemake.output[0], "w", encoding="utf8") as file:
            file.writelines(list_of_lines)    

    if trans == "LOW":
        list_of_lines[328] = "\n" + 'var_trans_pcap.LO(z,z_alias,"HVAC400KV")$(trans_links(z,z_alias,trans)) = trans_links_cap(z,z_alias,trans);'+"\n"
        with open(snakemake.output[0], "w", encoding="utf8") as file:
            file.writelines(list_of_lines)    


build_transmission(snakemake.wildcards.trans)