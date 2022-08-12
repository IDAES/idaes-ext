#include "read_data.h"
#include<string>
#include<iostream>
#include<fstream>

std::vector< std::vector<double> > read_data(comp_enum comp, test_data::data_set_enum data_set, double u_off, double h_off, double s_off){
  std::string comp_str = comp_enum_table[comp];
  std::string data_set_str, nd_str = "ND";
  std::ostringstream file_name_stream;
  std::istringstream row_stream;
  std::fstream filestream;
  char row_buf[1500];
  char col_buf[256];


  // Really should have used a map, oh well
  if(data_set == test_data::vapor_set){
    data_set_str = "vapor";
  }
  else if(data_set == test_data::liquid_set){
    data_set_str = "liquid";
  }
  else if(data_set == test_data::supercritical_set){
    data_set_str = "supercritical";
  }
  else if(data_set == test_data::saturated_set){
    data_set_str = "saturated";
  }
  else if(data_set == test_data::mixed_set){
    data_set_str = "mixed";
  }

  file_name_stream << "test_data/" << comp_str << "_" << data_set_str << "_data.csv";
  std::cout << "Reading: "<< file_name_stream.str() << std::endl;

  filestream.open(file_name_stream.str(), std::fstream::in);
  filestream.getline(row_buf, 1500); // ignore the header
  int nrows = 0;
  while (filestream.getline(row_buf, 1500)){
    ++nrows;
  }
  std::vector< std::vector<double> > rows(nrows, std::vector<double>(12));
  filestream.clear();
  filestream.seekg(0);
  int r = 0;
  int c = 0;
  filestream.getline(row_buf, 1500); // skip header
  while (filestream.getline(row_buf, 1500)){
    row_stream.str(row_buf);
    c = 0;
    while (row_stream.getline(col_buf, 256, ',')){
      if (!nd_str.compare(col_buf)){
        rows[r][c] = nan("no data");
      }
      else if (c == test_data::col_enum::h_col){
        rows[r][c] = std::stod(col_buf) - h_off;
      }
      else if (c == test_data::col_enum::u_col){
        rows[r][c] = std::stod(col_buf) - u_off;
      }
      else if (c == test_data::col_enum::s_col){
        rows[r][c] = std::stod(col_buf) - s_off;
      }
      else{
        rows[r][c] = std::stod(col_buf);
      }
      ++c;
      if(c == 12) break; // don't want phase string
    }
    ++r;
  }
  return rows;
}


void sort_sat( //read sat data, and assume the rows are: liq, vap, liq vap ... and even number of data
    comp_enum comp,
    test_data::data_set_enum data_set,
    std::vector< std::vector<double> > *liq_ptr,
    std::vector< std::vector<double> > *vap_ptr
){
  std::vector< std::vector<double> > all_data = read_data(comp, data_set);
  unsigned long int i, n = all_data.size()/2;  // beter be even
  liq_ptr->resize(n);
  vap_ptr->resize(n);
  for(i=0; i<n; ++i){
    liq_ptr->at(i) = all_data[2*i];
    vap_ptr->at(i) = all_data[2*i + 1];
  }
}
