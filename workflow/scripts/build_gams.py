"""
Change different parameters in the GAMS code
"""


def build_gams(year, varnewpcapQ, enable_fixed_ratios):
    """
    Crudely modify the GAMS code textfile, will break if line numbers change
    """
    with open(snakemake.input[0], "r", encoding="utf8") as file:
        list_of_lines = file.readlines()
    list_of_lines[
        45
    ] = f'$setglobal codefolderpath "{snakemake.params.sharedcodepath}"\n'
    list_of_lines[70] = f'$setglobal weather_yr "{year}"\n'
    list_of_lines[71] = f'$setglobal dem_yr "{year}"\n'
    #\n is not included below becasue line still continues
    list_of_lines[
        636
    ] = f"sum((z,h)$(hr2yr_map(yr,h)),demand(z,h))*{snakemake.params.co2intensity}"

    # Only add the constraints if enable_fixed_ratios is True:
    if enable_fixed_ratios:    
        list_of_lines.insert(443, 'eq_new_pcap_sub1' + "\n")
        list_of_lines.insert(444, 'eq_new_pcap_sub2' + "\n")
        list_of_lines.insert(445, 'eq_new_pcap_sub3' + "\n")
        list_of_lines.insert(446, 'eq_new_pcap_sub4' + "\n")
        list_of_lines.insert(571, f'eq_new_pcap_sub1 .. var_new_pcap("Solar") =E= {varnewpcapQ[0]} * sum(g, var_new_pcap(g));' + "\n")
        list_of_lines.insert(572, "\n")
        list_of_lines.insert(573, f'eq_new_pcap_sub2 .. var_new_pcap("Windonshore") =E= {varnewpcapQ[1]} * sum(g, var_new_pcap(g));' + "\n")
        list_of_lines.insert(574, "\n")
        list_of_lines.insert(575, f'eq_new_pcap_sub3 .. var_new_pcap("Windoffshore") =E= {varnewpcapQ[2]} * sum(g, var_new_pcap(g));' + "\n")
        list_of_lines.insert(576, "\n")
        list_of_lines.insert(577, f'eq_new_pcap_sub4 .. var_new_pcap("Windoffshorefloating") =E= {varnewpcapQ[3]} * sum(g, var_new_pcap(g));' + "\n")
        list_of_lines.insert(578, "\n")

    with open(snakemake.output[0], "w", encoding="utf8") as file:
        file.writelines(list_of_lines)


build_gams(snakemake.wildcards.year, snakemake.params.varnewpcapQ, snakemake.params.enable_fixed_ratios)
