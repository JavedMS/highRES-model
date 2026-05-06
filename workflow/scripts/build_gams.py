"""
Change different parameters in the GAMS code
"""
import json

def build_gams(year, varnewpcapQ, enable_fixed_ratios):
    """
    Crudely modify the GAMS code textfile, will break if line numbers change
    """
    with open(snakemake.input[0], "r", encoding="utf8") as file:
        list_of_lines = file.readlines()
    list_of_lines[
        46
    ] = f'$setglobal codefolderpath "{snakemake.params.sharedcodepath}"\n'
    list_of_lines[49] = '$setglobal gdx2sql "OFF" \n'
    list_of_lines[70] = f'$setglobal weather_yr "{year}"\n'
    list_of_lines[71] = f'$setglobal dem_yr "{year}"\n'
    #\n is not included below becasue line still continues
    list_of_lines[
        656
    ] = f"sum((z,h)$(hr2yr_map(yr,h)),demand(z,h))*{snakemake.params.co2intensity}"

    list_of_lines[541] = 'sum(gen_lim(z,g)$(not sameas(g, "Export")), var_gen(h,z,g))' + "\n"
    list_of_lines[560] = '=E= demand(z,h) + var_gen(h,z,"Export")$gen_lim(z,"Export");' + "\n"
    list_of_lines.insert(449, 'eq_imp_exp' + "\n")
    list_of_lines.insert(578,
    'eq_imp_exp .. sum((h,gen_lim(z,"Import")), var_gen(h,z, "Import")) =E= sum((h,gen_lim(z,"Export")), var_gen(h,z, "Export"));' + "\n")
    list_of_lines.insert(579, "\n")

    # Get the fylke_tech scenario data
    fylke_tech_key = snakemake.wildcards.fylke_tech
    fylke_tech_limit = json.loads(snakemake.params.fylke_tech_limit)
    scenario_data = fylke_tech_limit[fylke_tech_key]

    # Create parameter strings for solar and wind
    par_vre_new_pcap_s = ", ".join([f"{fylke} {values[0]}" for fylke, values in scenario_data.items()])
    par_vre_new_pcap_w = ", ".join([f"{fylke} {values[1]}" for fylke, values in scenario_data.items()])

    list_of_lines.insert(390, 'parameters' + "\n")
    list_of_lines.insert(391, f'par_vre_new_pcap_s(z) /{par_vre_new_pcap_s}/' + "\n")
    list_of_lines.insert(392, f'par_vre_new_pcap_w(z) /{par_vre_new_pcap_w}/;' + "\n")
    #list_of_lines.insert(385, 'par_vre_new_pcap_s(z) /NO03 1, NO11 1, NO15 1, NO18 1, NO30 1, NO34 1, NO38 1, NO42 1, NO46 1, NO50 1, NO54 1/' + "\n")
    #list_of_lines.insert(386, 'par_vre_new_pcap_w(z) /NO03 1, NO11 1, NO15 1, NO18 1, NO30 1, NO34 1, NO38 1, NO42 1, NO46 1, NO50 1, NO54 1/;' + "\n")
    list_of_lines.insert(453, 'eq_vre_new_pcap_s' + "\n")
    list_of_lines.insert(454, 'eq_vre_new_pcap_w' + "\n")
    list_of_lines.insert(583, 'eq_vre_new_pcap_s(z)$gen_lim(z,"Solar") .. var_new_pcap_z(z, "Solar") =L= par_vre_new_pcap_s(z)*var_new_pcap("Solar");'+"\n")
    list_of_lines.insert(584, "\n")
    #list_of_lines.insert(575, 'eq_vre_new_pcap_w(z)$gen_lim(z,"Windonshore") .. var_new_pcap_z(z, "Windonshore") =L= par_vre_new_pcap_w(z)*var_new_pcap("Windonshore");'+"\n")
    list_of_lines.insert(585,'eq_vre_new_pcap_w(z)$(gen_lim(z,"Windonshore_OF") or gen_lim(z,"Windonshore_F")) .. var_new_pcap_z(z, "Windonshore_OF") + var_new_pcap_z(z, "Windonshore_F") =L= par_vre_new_pcap_w(z)*(var_new_pcap("Windonshore_OF") + var_new_pcap("Windonshore_F"));'+"\n")
    list_of_lines.insert(586, "\n")

    #list_of_lines.insert(455, 'eq_store_PHS_limit' + "\n")
    #list_of_lines.insert(588, 'eq_store_PHS_limit .. sum((z,h)$s_lim(z,"PumpedHydro"),var_store_gen(h,z,"PumpedHydro")) =L= 1761;' + "\n")
    #list_of_lines.insert(589, "\n")
    #list_of_lines.insert(455, 'eq_new_pcap_sub1' + "\n")
    #list_of_lines.insert(588, 'eq_new_pcap_sub1 .. var_new_pcap("Solar") =L= 30;' + "\n")
    #list_of_lines.insert(589, "\n")
    # Only add the constraints if enable_fixed_ratios is True:
    if enable_fixed_ratios:    
        list_of_lines.insert(455, 'eq_new_pcap_sub1' + "\n")
        list_of_lines.insert(456, 'eq_new_pcap_sub2' + "\n")
        list_of_lines.insert(591, f'eq_new_pcap_sub1 .. var_new_pcap("Solar") =L= {varnewpcapQ[0]} * sum(g, var_new_pcap(g));' + "\n")
        list_of_lines.insert(592, "\n")
        list_of_lines.insert(593, f'eq_new_pcap_sub2 .. var_new_pcap("Windonshore_OF")+var_new_pcap("Windonshore_F")=L= {varnewpcapQ[1]} * sum(g, var_new_pcap(g));' + "\n")
        list_of_lines.insert(594, "\n")


    list_of_lines[159] = 'vre_lim(vre,z,r)=((area(vre,z,r)+(gen_exist_pcap_z(z,vre,"FX")$sameas(z,r)))>0.); \n'
    list_of_lines[165] = 'gen_lim(z,vre)=((sum(r,area(vre,z,r))+sum(lt,gen_exist_pcap_z(z,vre,lt)))>0.); \n'
    list_of_lines.insert(455, 'eq_costs_gen_vreconnection_bottom' + "\n")
    list_of_lines.insert(456, 'eq_costs_gen_vreconnection_floating' + "\n")

    list_of_lines.insert(273, 'costs_gen_vreconnection_bottom(z)' + "\n")
    list_of_lines.insert(274, 'costs_gen_vreconnection_floating(z)' + "\n")

    equation_text1 = (
    'eq_costs_gen_vreconnection_bottom(z) .. costs_gen_vreconnection_bottom(z) =E=\n'
    '    sum(vre_lim("Windoffshore",z,r), sum(z_alias$(sameas(z_alias,r) and trans_links_dist(z,z_alias,"HVACWindoffshore")),\n'
    '        var_new_vre_pcap_r(z,"Windoffshore",r)*((trans_links_dist(z,z_alias,"HVACWindoffshore")*trans_line_capex("HVACWindoffshore"))\n'
    '       +trans_sub_capex("HVACWindoffshore"))));' + "\n"
)

    equation_text2 = (
    'eq_costs_gen_vreconnection_floating(z) .. costs_gen_vreconnection_floating(z) =E=\n'
    '    sum(vre_lim("Windoffshorefloating",z,r), sum(z_alias$(sameas(z_alias,r) and trans_links_dist(z,z_alias,"HVDCWindoffshore")),\n'
    '        var_new_vre_pcap_r(z,"Windoffshorefloating",r)*((trans_links_dist(z,z_alias,"HVDCWindoffshore")*trans_line_capex("HVDCWindoffshore"))\n'
    '        +trans_sub_capex("HVDCWindoffshore"))));' + "\n"
)
    list_of_lines.insert(541, equation_text1)
    list_of_lines.insert(542, equation_text2)    

    list_of_lines.insert(492, '+costs_gen_vreconnection_bottom(z)' + "\n")
    list_of_lines.insert(493, '+costs_gen_vreconnection_floating(z)' + "\n")


    with open(snakemake.output[0], "w", encoding="utf8") as file:
        file.writelines(list_of_lines)


build_gams(snakemake.wildcards.year, snakemake.params.varnewpcapQ, snakemake.params.enable_fixed_ratios)
#build_gams(snakemake.wildcards.year)
