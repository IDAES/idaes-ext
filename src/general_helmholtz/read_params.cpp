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
    std::ostringstream nl_file_path, json_file_path, except_string;

    if(data_path == ""){
        if(const char* env_data_path = getenv("IDAES_HELMHOLTZ_DATA_PATH")){
            data_path = env_data_path;
        }
    }
    
    // IDAES component naming standards are not clear on case.  I'll assume you
    // won't have two different components whoes name differs only by case. Some
    // file systems aren't case sensitive anyway.  This setup requires the file
    // names for the components to be lower case, but the compoent name can be
    // upper, lower, or mixed case.  Like I could have "h2o" and/or "H2O", but
    // both componenets would use the "h2o" parameter set.

    //Get the lower case component name for file names.
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

    for(uint i=0; i<14; i++){
        cdata[comp_idx].expr_map[i] = boost::json::value_to<long>(jp.at("expr_map").at(i));
    }
    for(uint i=0; i<3; i++){
        cdata[comp_idx].var_map[i] = boost::json::value_to<long>(jp.at("var_map").at(i));
    }

    nl_file_path << data_path << boost::json::value_to<std::string>(jp.at("nl_file"));
    std::string nl_file_string = nl_file_path.str();
    std::ifstream nl_file_stream(nl_file_path.str());
    if(!nl_file_stream.good()){
        except_string << "Data file not found: " << nl_file_path.str();
        std::cout << except_string.str() << std::endl;
        return MISSING_DATA;
    }
    // Read the NL-file
    asl = ASL_alloc(ASL_read_pfgh);
    pfgh_read(jac0dim(nl_file_string.c_str(), nl_file_string.length()), 0);
    cdata[comp_idx].asl = (void*)asl;
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
    if (it->value().as_object().contains("u_off")){
        test_data.at(i).u_off = boost::json::value_to<double>(it->value().at("u_off"));
    }
    if (it->value().as_object().contains("h_off")){
        test_data.at(i).h_off = boost::json::value_to<double>(it->value().at("h_off"));
    }
    if (it->value().as_object().contains("s_off")){
        test_data.at(i).s_off = boost::json::value_to<double>(it->value().at("s_off"));
    }
    ++i;
  }

  return test_data;
}