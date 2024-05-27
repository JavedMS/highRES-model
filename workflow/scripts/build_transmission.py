"""
Change transmission node capacities in the GAMS code
"""


def build_transmission(trans):
    """
    Crudely modify the GAMS code textfile, will break if line numbers change
    """

    if trans == "DefaultYes":
        with open(snakemake.input[0], "r", encoding="utf8") as file:
            list_of_lines = file.readlines()
        list_of_lines[72] = '$setglobal fx_trans "YES"'+"\n"        

        with open(snakemake.output[0], "w", encoding="utf8") as file:
            file.writelines(list_of_lines)


    if trans == "DefaultNO":
        with open(snakemake.input[0], "r", encoding="utf8") as file:
            list_of_lines = file.readlines()
        list_of_lines[72] = '$setglobal fx_trans "NO"'+"\n"

        with open(snakemake.output[0], "w", encoding="utf8") as file:
            file.writelines(list_of_lines)

    if trans == "FIXEDOH":
        with open(snakemake.input[0], "r", encoding="utf8") as file:
            list_of_lines = file.readlines()
        list_of_lines[72] = '$setglobal fx_trans "NO"'+"\n"        
        list_of_lines[328] = "\n" + 'var_trans_pcap.FX(z,z_alias,"HVAC400KV")$(trans_links(z,z_alias,"HVAC400KV")) = trans_links_cap(z,z_alias,"HVAC400KV");'+"\n"
        list_of_lines[329] = 'var_trans_pcap.UP(z,z_alias,"HVDCSubsurface")$(trans_links(z,z_alias,"HVDCSubsurface")) = 50.;'+"\n"
        
        with open(snakemake.output[0], "w", encoding="utf8") as file:
            file.writelines(list_of_lines)    

    if trans == "LOWOH":
        with open(snakemake.input[0], "r", encoding="utf8") as file:
            list_of_lines = file.readlines()
        list_of_lines[72] = '$setglobal fx_trans "NO"'+"\n"        
        list_of_lines[328] = "\n" + 'var_trans_pcap.LO(z,z_alias,"HVAC400KV")$(trans_links(z,z_alias,"HVAC400KV")) = trans_links_cap(z,z_alias,"HVAC400KV");'+"\n"
        with open(snakemake.output[0], "w", encoding="utf8") as file:
            file.writelines(list_of_lines)    


build_transmission(snakemake.wildcards.trans)