#include<boost/json/src.hpp>
#include<cstdlib>
#include<math.h>
#include"config.h"
#include"read_params.h"
#include<iostream>
#include<fstream>
#include<string>
#include<asl.h>

std::unordered_map<std::string, uint> cindex;
std::vector<parameters_struct> cdata;

uint read_params(std::string comp, std::string data_path){
    try{
        return cindex.at(comp);
    }
    catch(std::out_of_range const&){
    } // not read yet
    ASL *asl;
    std::ostringstream 
        nl_file_path, 
        nl_file_tcx_path, 
        nl_file_visc_path, 
        nl_file_st_path, 
        json_file_path, 
        except_string;
    if(data_path == ""){
        if(const char* env_data_path = getenv("IDAES_HELMHOLTZ_DATA_PATH")){
            data_path = env_data_path;
        }
    } else{
        // For now at least I have to assume all the data files are in the same place, since
        // some expressions defined in NL files may call functions in the property library
        // and this may spawn another process, I need to transfer the data path set in the
        // function call to the environment variable to ensure that the data files can be
        // located by the new process.
        setenv("IDAES_HELMHOLTZ_DATA_PATH", data_path.c_str(), 1);
    }
    // I'll assume component names are not case sensitive, so get lower case name.
    std::string lower_comp = comp;
    for(uint i = 0; i < comp.length(); i++){
        lower_comp[i] = tolower(comp[i]);
    }
    // parameter file other parameters
    json_file_path << data_path << lower_comp << "_parameters.json";
    // Read the json parameter file
    std::ifstream jp_file_stream(json_file_path.str());
    if(!jp_file_stream.good()){
        except_string << "Data file not found: " << json_file_path.str();
        std::cout << except_string.str() << std::endl;
        return MISSING_DATA;
    }
    std::stringstream jp_buffer;
    jp_buffer << jp_file_stream.rdbuf();
    boost::json::value jp = boost::json::parse(jp_buffer.str());

    parameters_struct _empty;
    cindex[comp] = cdata.size();
    cdata.push_back(_empty);
    uint comp_idx = cindex[comp];

    cdata[comp_idx].MW = boost::json::value_to<double>(jp.at("param").at("MW"));
    cdata[comp_idx].R = boost::json::value_to<double>(jp.at("param").at("R"));
    cdata[comp_idx].T_star = boost::json::value_to<double>(jp.at("param").at("T_star"));
    cdata[comp_idx].rho_star = boost::json::value_to<double>(jp.at("param").at("rho_star"));
    cdata[comp_idx].Tc = boost::json::value_to<double>(jp.at("param").at("Tc"));
    cdata[comp_idx].rhoc = boost::json::value_to<double>(jp.at("param").at("rhoc"));
    cdata[comp_idx].Pc = boost::json::value_to<double>(jp.at("param").at("Pc"));
    cdata[comp_idx].Tt = boost::json::value_to<double>(jp.at("param").at("Tt"));
    cdata[comp_idx].Pt = boost::json::value_to<double>(jp.at("param").at("Pt"));
    cdata[comp_idx].rhot_l = boost::json::value_to<double>(jp.at("param").at("rhot_l"));
    cdata[comp_idx].rhot_v = boost::json::value_to<double>(jp.at("param").at("rhot_v"));
    cdata[comp_idx].P_min = boost::json::value_to<double>(jp.at("param").at("P_min"));
    cdata[comp_idx].P_max = boost::json::value_to<double>(jp.at("param").at("P_max"));
    cdata[comp_idx].rho_max = boost::json::value_to<double>(jp.at("param").at("rho_max"));
    cdata[comp_idx].T_min = boost::json::value_to<double>(jp.at("param").at("T_min"));
    cdata[comp_idx].T_max = boost::json::value_to<double>(jp.at("param").at("T_max"));
    try{
        cdata[comp_idx].rs_offset1 = boost::json::value_to<double>(jp.at("param").at("reference_state_offset").at(0));
        cdata[comp_idx].rs_offset2 = boost::json::value_to<double>(jp.at("param").at("reference_state_offset").at(1));
    }
    catch(std::out_of_range const&){
        cdata[comp_idx].rs_offset1 = 0;
        cdata[comp_idx].rs_offset2 = 0;
    }
    for(uint i=0; i<expr_map_size; i++){
        cdata[comp_idx].expr_map[i] = boost::json::value_to<long>(jp.at("expr_map").at(i));
    }
    for(uint i=0; i<var_map_size; i++){
        cdata[comp_idx].var_map[i] = boost::json::value_to<long>(jp.at("var_map").at(i));
    }
    cdata[comp_idx].have_tcx = boost::json::value_to<bool>(jp.at("have_tcx"));
    cdata[comp_idx].have_visc = boost::json::value_to<bool>(jp.at("have_visc"));
    cdata[comp_idx].have_st = boost::json::value_to<bool>(jp.at("have_st"));
    nl_file_path << data_path << boost::json::value_to<std::string>(jp.at("nl_file"));
    std::string nl_file_string = nl_file_path.str();
    // Read the NL-file
    asl = ASL_alloc(ASL_read_pfgh);
    pfgh_read(jac0dim(nl_file_string.c_str(), nl_file_string.length()), 0);
    cdata[comp_idx].asl = (void*)asl;

    if(cdata[comp_idx].have_tcx){
        for(uint i=0; i<var_map_size; i++){
            cdata[comp_idx].var_map_tcx[i] = boost::json::value_to<long>(jp.at("var_map_tcx").at(i));
        }
        nl_file_path.clear();
        nl_file_tcx_path << data_path << boost::json::value_to<std::string>(jp.at("nl_file_tcx"));
        nl_file_string = nl_file_tcx_path.str();
        asl = ASL_alloc(ASL_read_pfgh);
        pfgh_read(jac0dim(nl_file_string.c_str(), nl_file_string.length()), 0);
        cdata[comp_idx].asl_tcx = (void*)asl;
    }
    if(cdata[comp_idx].have_visc){
        for(uint i=0; i<var_map_size; i++){
            cdata[comp_idx].var_map_visc[i] = boost::json::value_to<long>(jp.at("var_map_visc").at(i));
        }
        nl_file_path.clear();
        nl_file_visc_path << data_path << boost::json::value_to<std::string>(jp.at("nl_file_visc"));
        nl_file_string = nl_file_visc_path.str();
        asl = ASL_alloc(ASL_read_pfgh);
        pfgh_read(jac0dim(nl_file_string.c_str(), nl_file_string.length()), 0);
        cdata[comp_idx].asl_visc = (void*)asl;
    }
    if(cdata[comp_idx].have_st){
        for(uint i=0; i<var_map_size; i++){
            cdata[comp_idx].var_map_st[i] = boost::json::value_to<long>(jp.at("var_map_st").at(i));
        }
        nl_file_path.clear();
        nl_file_st_path << data_path << boost::json::value_to<std::string>(jp.at("nl_file_st"));
        nl_file_string = nl_file_st_path.str();
        asl = ASL_alloc(ASL_read_pfgh);
        pfgh_read(jac0dim(nl_file_string.c_str(), nl_file_string.length()), 0);
        cdata[comp_idx].asl_st = (void*)asl;
    }
    return comp_idx;
}

std::vector<tests_struct> read_run_tests(void){
  std::vector<tests_struct> test_data;
  std::ostringstream file_name_stream, except_string;
  std::string data_path("");
  
  if(const char* env_data_path = getenv("IDAES_HELMHOLTZ_TEST_DATA_PATH")){
    data_path = env_data_path;
  }
  
  file_name_stream << data_path << "run_tests.json";
  std::cout << "Reading: "<< file_name_stream.str() << std::endl;
  std::ifstream jp_file_stream(file_name_stream.str());
  
  if(!jp_file_stream.good()){
    except_string << "Data file not found: " << file_name_stream.str();
    std::cout << except_string.str() << std::endl;
    throw std::runtime_error(except_string.str());
  }

  std::stringstream jp_buffer;
  jp_buffer << jp_file_stream.rdbuf();

  boost::json::value jp = boost::json::parse(jp_buffer.str());

  test_data.resize(jp.as_object().size());
  uint i = 0;
  for (auto it=jp.as_object().cbegin(); it != jp.as_object().cend(); ++it){
    std::cout << i << ": " << it->key() << std::endl;
    test_data.at(i).comp_str = boost::json::value_to<std::string>(it->value().at("component"));
    test_data.at(i).test_set = boost::json::value_to<std::string>(it->value().at("tests"));
    ++i;
  }
  return test_data;
}

void set_reference_state_offset(uint comp_idx, double n1, double n2){
    cdata[comp_idx].rs_offset1 = n1;
    cdata[comp_idx].rs_offset2 = n2;
}

