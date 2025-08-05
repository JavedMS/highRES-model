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
        list_of_lines[73] = '$setglobal trans_inv "OFF"'+"\n"        

        with open(snakemake.output[0], "w", encoding="utf8") as file:
            file.writelines(list_of_lines)


    if trans == "LOWOH":
        with open(snakemake.input[0], "r", encoding="utf8") as file:
            list_of_lines = file.readlines()
        list_of_lines[73] = '$setglobal trans_inv "USER"'+"\n"
        list_of_lines[74] = '$setglobal trans_cap_lim "20"'+"\n"        #cap value to be set by user 

        with open(snakemake.output[0], "w", encoding="utf8") as file:
            file.writelines(list_of_lines)

    if trans == "FIXEDOH":
        with open(snakemake.input[0], "r", encoding="utf8") as file:
            list_of_lines = file.readlines()
        list_of_lines[73] = '$setglobal trans_inv "USER"'+"\n"
        list_of_lines[74] = '$setglobal trans_cap_lim "20"'+"\n"        #cap value to be set by user        
        list_of_lines[334] = 'var_new_trans_pcap.UP(z,z_alias,"HVAC400KV")$(trans_links(z,z_alias,"HVAC400KV")) = 0;'+"\n"
        list_of_lines[335] = 'var_new_trans_pcap.UP(z,z_alias,"HVDCSubsurface")$(trans_links(z,z_alias,"HVDCSubsurface")) = %trans_cap_lim%;'+"\n"
        
        with open(snakemake.output[0], "w", encoding="utf8") as file:
            file.writelines(list_of_lines)    

    if trans == "SubsurfOFF":
        with open(snakemake.input[0], "r", encoding="utf8") as file:
            list_of_lines = file.readlines()
        list_of_lines[73] = '$setglobal trans_inv "USER"'+"\n"
        list_of_lines[74] = '$setglobal trans_cap_lim "20"'+"\n"        #cap value to be set by user        
        list_of_lines[334] = 'var_new_trans_pcap.UP(z,z_alias,"HVAC400KV")$(trans_links(z,z_alias,"HVAC400KV")) = %trans_cap_lim%;'+"\n"
        list_of_lines[335] = 'var_new_trans_pcap.UP(z,z_alias,"HVDCSubsurface")$(trans_links(z,z_alias,"HVDCSubsurface")) = 0;'+"\n"
        with open(snakemake.output[0], "w", encoding="utf8") as file:
            file.writelines(list_of_lines)    

    if trans == "Subsurfmust":
        with open(snakemake.input[0], "r", encoding="utf8") as file:
            list_of_lines = file.readlines()
        list_of_lines[72] = '$setglobal fx_trans "NO"'+"\n" 
        list_of_lines[328] = 'var_trans_pcap.LO(z,z_alias,"HVAC400KV")$(trans_links(z,z_alias,"HVAC400KV")) = trans_links_cap(z,z_alias,"HVAC400KV");'+"\n"       
        list_of_lines.insert(451, 'eq_new_trans_pcap' + "\n")
        #list_of_lines.insert(637, "\n")
        list_of_lines.insert(637,
        f'eq_new_trans_pcap .. sum(trans_links(z,z_alias,"HVDCSubsurface"), var_trans_pcap(z,z_alias,"HVDCSubsurface")) =G= '
        f'0.55 * (sum(trans_links(z,z_alias,trans), var_trans_pcap(z,z_alias,trans))-sum(trans_links(z,z_alias,trans), trans_links_cap(z,z_alias,trans)));' + "\n")
        list_of_lines.insert(638, "\n")
        with open(snakemake.output[0], "w", encoding="utf8") as file:
            file.writelines(list_of_lines)


build_transmission(snakemake.wildcards.trans)