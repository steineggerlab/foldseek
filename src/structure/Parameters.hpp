#pragma once
#include <iostream>
#include <cstring>
#include <algorithm>
#include <cctype>
#include <filesystem>
#include <vector>

using namespace std;
namespace fs = std::filesystem;

class Parameters{
    private:
        bool show_structure = false;
        bool no_panel = false;
        bool arg_okay = true;
        bool benchmark_mode = false;
        vector<string> in_file;
        string chainfile = "";
        string mode = "protein";
        string msa_file = "";
        string foldmason_file = "";
        // -fst : target 구조 소스 (Foldseek DB | 구조 디렉터리 | 구조 파일 | "auto")
        string foldseek_target = "";
        // -fsr : Foldseek 결과 (m8 12/17/21/29 컬럼 | 멀티머 _report 14 컬럼)
        string foldseek_result = "";
    public:
        Parameters(int argc, char* argv[]);

        void print_args();

        // get, set
        vector<string>& get_in_file(){
            return in_file;
        }
        string get_in_file(int idx){
            if (idx < in_file.size()){
                return in_file[idx];
            }
            return "";
        }
        string get_chainfile(){
            return chainfile;
        }
        string get_mode(){
            return mode;
        }
        bool get_show_structure(){
            return show_structure;
        }
        bool get_no_panel(){
            return no_panel;
        }
        bool get_benchmark_mode(){
            return benchmark_mode;
        }
        bool check_arg_okay(){
            return arg_okay;
        }
        string get_msa_file(){
            return msa_file;
        }
        string get_foldmason_file(){
            return foldmason_file;
        }
        string get_foldseek_target(){
            return foldseek_target;
        }
        string get_foldseek_result(){
            return foldseek_result;
        }
};