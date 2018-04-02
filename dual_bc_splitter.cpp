#include <string>
#include <vector>
#include <fstream>
#include <iostream>
#include <cstdlib> 
#include <cstring>
#include <cctype>
#include <sstream>
#include <algorithm>
#include <memory>
#include <tuple>
#include <map>
#include <vector>
#include <ctime>
#include <iomanip>
#include <cassert>

#include "BKTree.h"
#include "fastq_reader.hpp"
#include "fastq_writer.hpp"
#include "barcode_loader.hpp"

#include <boost/regex.hpp>
#include <boost/program_options.hpp>
#include <boost/algorithm/string.hpp>

#include <sys/stat.h>
#include <sys/types.h>
#include <unistd.h>

namespace po = boost::program_options;

typedef std::tuple<std::string, std::string, int, std::string> match_t;

class dual_bc_splitter {
    public:
        bool parse_args(int argc, char* argv[]);
        void initialize();
        void core_engine();
        void print_help();
        void load_barcodes();
        std::string vec_to_str(std::vector<std::string> my_vector);
        void write_outfile(std::string outfile);
        std::string get_key_val(unsigned int bc1_pos, unsigned int bc2_pos);

        unsigned int updateMaps(std::string barcode_str,
            std::vector<std::string>& bc1_words, std::vector<std::string>& bc2_words,
            std::vector<std::string>& r1_words, std::vector<std::string>& r2_words);
        void write_for_single_file(std::string& barcode,
          std::map<std::string, std::vector<std::unique_ptr<std::string>>>& readMap, std::ofstream& lof);
        void get_lines_to_vec(fastq_reader& lreader, std::vector<std::string>& lvec, 
            int lcount);
        std::string write_entry(const match_t& bc1_tuple, const match_t& bc2_tuple);
        void writeMapsToFile(std::ofstream& count_log);
        bool has_suffix(const std::string &str, const std::string &suffix);
        void compress_files();
    private:
		// Input variables pushed through command line
		std::string bc1_file_str;
		std::string bc2_file_str;
        std::string bc1_read_str;
        std::string bc2_read_str;
		std::string read1_str;
		std::string read2_str;
		std::string prefix_str;
		std::string outdirpath;
		int mm_bc1;
		int mm_bc2;
        long total_run = 0;
        unsigned int allowed_MB;
  
        unsigned long unique_exact_count = 0;
        unsigned long unique_nonexact_count = 0;
        unsigned long ambiguous_count = 0;
        unsigned long no_match_count = 0;
        unsigned long total_count = 0;
        

        po::options_description desc;

        BCLoader bc1_loader;
        BCLoader bc2_loader;
        std::unordered_map<std::string, unsigned int> count_map_exact; 
        std::unordered_map<std::string, unsigned int> count_map_nonexact; 
        std::set<std::string> outfile_set;
        std::set<std::string> outfile_all_set;

        std::map<std::string, std::vector<std::unique_ptr<std::string>>> r1QMap;
        std::map<std::string, std::vector<std::unique_ptr<std::string>>> r2QMap;
        std::map<std::string, std::vector<std::unique_ptr<std::string>>> bc1QMap;
        std::map<std::string, std::vector<std::unique_ptr<std::string>>> bc2QMap;

        std::map<std::string, std::unique_ptr<fastq_writer>> read1_writer_map;
        std::map<std::string, std::unique_ptr<fastq_writer>> read2_writer_map;
        std::map<std::string, std::unique_ptr<fastq_writer>> barcode1_writer_map;
        std::map<std::string, std::unique_ptr<fastq_writer>> barcode2_writer_map;

};


void dual_bc_splitter::print_help() {
    std::cout << desc << "\n";
	std::cout << "Usage: concensus --bc1 <bc1 file> --bc2 <bc2 file> "
        "--file1 <file1> --file2 <file2> "
        " -p <prefix_str> -o <outdir>\n\n";
}

std::string dual_bc_splitter::get_key_val(unsigned int bc1_pos, unsigned int bc2_pos) {
    std::string bc1_pos_str = std::to_string(bc1_pos);
    std::string bc2_pos_str = std::to_string(bc2_pos);
    std::string combined = bc1_pos_str + "_" + bc2_pos_str;
    return combined;
}

void dual_bc_splitter::load_barcodes() {

    bc1_loader = BCLoader(bc1_file_str);
    bc1_loader.load_map();
    bc1_loader.load_tree();
    bc1_loader.print_map();
    //bc1_loader.print_name_to_index();
    //bc1_loader.print_seq_to_index();

    bc2_loader = BCLoader(bc2_file_str);
    bc2_loader.load_map();
    bc2_loader.load_tree();
    bc2_loader.print_map();
    //bc2_loader.print_name_to_index();
    //bc2_loader.print_seq_to_index();

}

bool dual_bc_splitter::parse_args(int argc, char* argv[]) {
    
    bool all_set = true;
    

    desc.add_options()
		("help,h", "produce help message")
		("bc1", po::value<std::string>(&bc1_file_str), "barcode 1 dict file")
		("bc2", po::value<std::string>(&bc2_file_str), "barcode 2 dict file")
		("bc1_read", po::value<std::string>(&bc1_read_str), "barcode 1 fastq file")
		("bc2_read", po::value<std::string>(&bc2_read_str), "barcode 2 fastq file")
		("read1", po::value<std::string>(&read1_str), "read 1 fastq file")
		("read2", po::value<std::string>(&read2_str), "read 2 fastq file")
		("prefix,p", po::value<std::string>(&prefix_str), "Prefix string")
		("outdir,o", po::value<std::string>(&outdirpath), "Output directory")	
		("mm_bc1", po::value(&mm_bc1)->default_value(1), "Optional/Maximum allowed mismatches for barcode 1.")
		("mm_bc2", po::value(&mm_bc2)->default_value(1), "Optional/Maximum allowed mismatches for barcode 2.")
        ("total_run", po::value(&total_run)->default_value(-1), "Optional/total reads to process.")
        ("allowed-mb", po::value(&allowed_MB)->default_value(2048),
            "Optional/Estimated memory requirement in MB.")
	;

    po::variables_map vm;
    po::store(po::parse_command_line(argc, argv, desc), vm);
    po::notify(vm);

    if (vm.count("help")) {
        //print_help();
        return 0;
    } else {
        //all_set = false;
    }

    if (vm.count("bc1")) {
		std::cout << "barcode 1 dict file is set to: " << bc1_file_str << ".\n";
	} else {
		all_set = false;
		std::cout << "Error: barcode 1 dict file is not set.\n";
	}

    if (vm.count("bc2")) {
		std::cout << "barcode 2 dict file is set to: " << bc2_file_str << ".\n";
	} else {
		all_set = false;
		std::cout << "Error: barcode 2 dict file is not set.\n";
	}

	if (vm.count("bc1_read")) {
		std::cout << "bc1 fastq file is set to: " << bc1_read_str << ".\n";
	} else {
		all_set = false;
		std::cout << "Error: bc1 fastq file is not set.\n";
	}

    if (vm.count("bc2_read")) {
		std::cout << "bc2 fastq file is set to: " << bc2_read_str << ".\n";
	} else {
		all_set = false;
		std::cout << "Error: bc2 fastq file is not set.\n";
	}

    if (vm.count("read1")) {
		std::cout << "Read 1 fastq file is set to: " << read1_str << ".\n";
	} else {
		all_set = false;
		std::cout << "Error: read 1 fastq file is not set.\n";
	}

	if (vm.count("read2")) {
		std::cout << "Read 2 fastq file is set to: " << read2_str << ".\n";
	} else {
		all_set = false;
		std::cout << "Error: read 2 fastq file is not set.\n";
	}    


    if (vm.count("prefix")) {
		std::cout << "Prefix string is set to: " << prefix_str << ".\n";
	} else {
		std::cout << "Error: Prefix string is not set.\n";
	}
 
    if (vm.count("outdir")) {
		std::cout << "Outdir is set to: " << outdirpath << ".\n";
	} else {
		all_set = false;
		std::cout << "Error: Outdir is not set.\n";
	}

    std::cout << "Maximum allowed mismatches for barcode 1 is set to " << mm_bc1 << ".\n";
    std::cout << "Maximum allowed mismatches for barcode 2 is set to " << mm_bc2 << ".\n";
    return all_set;

}

void dual_bc_splitter::initialize() {
        
	struct stat st = {0};

	if (stat(outdirpath.c_str(), &st) == -1) {
		mkdir(outdirpath.c_str(), 0755);
	}
}


std::string dual_bc_splitter::vec_to_str(std::vector<std::string> lvec) {
    std::string res_str = "";
    for (int j = 0; j < lvec.size(); j++) {
        if (j == 0) {
            res_str += lvec[j];
        } else {
            res_str += ", " + lvec[j];
        }
    }
    return res_str;
    
}

std::string dual_bc_splitter::write_entry(const match_t& bc1_tuple, const match_t& bc2_tuple) {

    std::string no_match_str = "no_match";
    std::string ambiguous_str = "ambiguous";
    std::string unique_str = "unique";
    std::string retval;

    std::string bc1_mt = std::get<0>(bc1_tuple);
    std::string bc2_mt = std::get<0>(bc2_tuple);

    if (0 == bc1_mt.compare(no_match_str) || 
            0 == bc2_mt.compare(no_match_str)) {

        no_match_count++;
        retval = "no_match";

    } else if (0 == bc1_mt.compare(ambiguous_str) || 
            0 == bc2_mt.compare(ambiguous_str)) {

        ambiguous_count++;
        retval = "ambiguous";
        
    } else if (0 == bc1_mt.compare(unique_str) && 
               0 == bc2_mt.compare(unique_str)) {

        std::string bc1_val = std::get<1>(bc1_tuple);
        std::string bc2_val = std::get<1>(bc2_tuple);

        retval = bc1_val + "_" + bc2_val;
        unsigned int bc1_index = bc1_loader.get_seq_to_index(bc1_val);
        unsigned int bc2_index = bc2_loader.get_seq_to_index(bc2_val);

        // None of the index should be zero
        assert(bc1_index > 0);
        assert(bc2_index > 0);

        std::string key_val = get_key_val(bc1_index, bc2_index);

        //std::cout << "All unique" << "\n";
        int bc1_dist = std::get<2>(bc1_tuple);         
        int bc2_dist = std::get<2>(bc2_tuple);         
 
        if (0 == bc1_dist && 
            0 == bc2_dist) {

            unique_exact_count++;
            count_map_exact[key_val] += 1;

        } else {   

            unique_nonexact_count++;  
            count_map_nonexact[key_val] += 1;

        }

    } else {
        std::string mystr = "Illegal match type(s): bc1_mt: " + bc1_mt + 
            " bc2_mt: " + bc2_mt;
        throw my_exception(mystr);
    }
    return retval;
        
}

unsigned int dual_bc_splitter::updateMaps(std::string barcode_str, 
    std::vector<std::string>& bc1_words, std::vector<std::string>& bc2_words,
    std::vector<std::string>& r1_words, std::vector<std::string>& r2_words) {

    unsigned int allcap = 0;
    for (auto lstr: bc1_words) {
        bc1QMap[barcode_str].push_back(std::make_unique<std::string>(lstr));
        allcap += lstr.capacity();
    }
    for (auto lstr: bc2_words) {
        bc2QMap[barcode_str].push_back(std::make_unique<std::string>(lstr));
        allcap += lstr.capacity();
    }
    for (auto lstr: r1_words) {
        r1QMap[barcode_str].push_back(std::make_unique<std::string>(lstr));
        allcap += lstr.capacity();
    }
    for (auto lstr: r2_words) {
        r2QMap[barcode_str].push_back(std::make_unique<std::string>(lstr));
        allcap += lstr.capacity();
    }


    return allcap;
}

void dual_bc_splitter::write_outfile(std::string outfile) {
    std::ofstream outf(outfile);
    if (outf.is_open()) {
        // Write header
        outf << "bc1_index" << "\t" << "bc2_index" << "\t" << "bc1_bc2_seq" << "\t" 
             << "read_count_exact" << "\t" << "read_count_nonexact" 
             << "\t" << "read_count_all" << "\n";
        std::vector<std::string> bc2_vec = bc2_loader.get_name_vector();
        std::vector<std::string> bc1_vec = bc1_loader.get_name_vector();

        for (const auto& bc1_ : bc1_vec) {
            for (const auto& bc2_ : bc2_vec) {
                unsigned int bc1_index = bc1_loader.get_name_to_index(bc1_);
                unsigned int bc2_index = bc2_loader.get_name_to_index(bc2_);

                std::string key_val = get_key_val(bc1_index, bc2_index);

                unsigned int val_exact = count_map_exact[key_val];
                unsigned int val_nonexact = count_map_nonexact[key_val];
                unsigned int val_all = val_exact + val_nonexact;
                std::string val_exact_str = std::to_string(val_exact);
                std::string val_nonexact_str = std::to_string(val_nonexact);
                std::string val_all_str = std::to_string(val_all);
                std::string bc1_val = bc1_loader.val_from_bc_map(bc1_);
                std::string bc2_val = bc2_loader.val_from_bc_map(bc2_);
                std::string bc1_bc2_seq = bc1_val + "_" + bc2_val;
                std::string out_str = bc1_ + "\t" + bc2_ + "\t" + bc1_bc2_seq +
                    "\t" + val_exact_str + "\t" + val_nonexact_str + "\t" + 
                    val_all_str;
                outf << out_str << "\n";   
            }
        } 
        outf << "ambiguous\t" << ambiguous_count << "\n";
        outf << "no_match\t" << no_match_count << "\n";

        outf.close();
    } else {
        throw my_exception("Problem in opening the output file.\n");
    }

    
}

void dual_bc_splitter::write_for_single_file(std::string& barcode, 
    std::map<std::string, std::vector<std::unique_ptr<std::string>>>& readMap,
    std::ofstream& lof) {
    
    std::vector<std::unique_ptr<std::string>> & valSet = readMap[barcode];
    for (auto const& kv : valSet) {
        std::string val = *kv;
        lof << val << "\n";
    }
}

void dual_bc_splitter::writeMapsToFile(std::ofstream& count_log) {

    for (auto& kv : r1QMap) {
        std::string barcode = kv.first;
        std::string file1;
        std::string file2;
        std::string bcfile1;
        std::string bcfile2;

        std::string unified_barcode;
        if (barcode.compare("no_match") == 0 ||
            barcode.compare("ambiguous") == 0 ) {
            unified_barcode = "unmatched";
          
        } else {
            unified_barcode = barcode;
        }

        if (unified_barcode.compare("unmatched") == 0) {
            file1 = outdirpath + "/" + prefix_str + ".unmatched.1.fastq";
            file2 = outdirpath + "/" + prefix_str + ".unmatched.2.fastq";
            bcfile1 = outdirpath + "/" + prefix_str + ".unmatched.barcode_1.fastq";
            bcfile2 = outdirpath + "/" + prefix_str + ".unmatched.barcode_2.fastq";
        
        } else {
            file1 = outdirpath + "/" + prefix_str + "." + barcode + ".unmapped.1.fastq";
            file2 = outdirpath + "/" + prefix_str + "." + barcode + ".unmapped.2.fastq";
            bcfile1 = outdirpath + "/" + prefix_str + "." + barcode + ".unmapped.barcode_1.fastq";
            bcfile2 = outdirpath + "/" + prefix_str + "." + barcode + ".unmapped.barcode_2.fastq";
        }

        std::ios_base::openmode lmode;
        if (outfile_set.count(unified_barcode) == 0) {
            outfile_set.insert(unified_barcode);
            lmode = std::ofstream::out|std::ofstream::trunc;
        } else {
            lmode = std::ofstream::out|std::ofstream::app;
        }

        // Note that we are basically writing simple text file, not gzipped
        // file, since software for efficiently appending to a gz file is 
        // pretty much not available.

        std::ofstream of_r1(file1, lmode);
        std::ofstream of_r2 (file2, lmode);
        std::ofstream of_bc1(bcfile1, lmode);
        std::ofstream of_bc2(bcfile2, lmode);

        outfile_all_set.insert(file1);
        outfile_all_set.insert(file2);
        outfile_all_set.insert(bcfile1);
        outfile_all_set.insert(bcfile2);

        write_for_single_file(barcode, r1QMap, of_r1);
        write_for_single_file(barcode, r2QMap, of_r2);
        write_for_single_file(barcode, bc1QMap, of_bc1);
        write_for_single_file(barcode, bc2QMap, of_bc2);
        count_log << "barcode: " << barcode << " bc1_written_reads: " << bc1QMap[barcode].size()/4 << "\n";
        count_log << "barcode: " << barcode << " bc2_written_reads: " << bc2QMap[barcode].size()/4 << "\n";
        count_log << "barcode: " << barcode << " r1_written_reads: " << r1QMap[barcode].size()/4 << "\n";
        count_log << "barcode: " << barcode << " r2_written_reads: " << r2QMap[barcode].size()/4 << "\n";
    }

    // Clear all contents of the read maps
    r1QMap.clear();
    r2QMap.clear();
    bc1QMap.clear();
    bc2QMap.clear();
}

void dual_bc_splitter::get_lines_to_vec(fastq_reader& lreader, 
    std::vector<std::string>& lvec, int lcount) {
    
    for (int j = 0; j < lcount; j++) {
            std::string lword;
            if (!lreader.getline(lword)) {
                std::string err_msg = "Truncated file: " + lreader.get_filename();
                throw my_exception(err_msg);
            }
            lvec.push_back(std::move(lword));
        }
}
void dual_bc_splitter::core_engine() {


    fastq_reader bc1_reader(bc1_read_str);
    fastq_reader bc2_reader(bc2_read_str);
    fastq_reader reader1(read1_str);
    fastq_reader reader2(read2_str);


    std::string count_log_file = outdirpath + "/" + prefix_str + "_logfile.txt";
    std::ofstream count_log(count_log_file);

    const unsigned long MB_SIZE = 1024 * 1024;
    unsigned long total_allowed = allowed_MB * MB_SIZE;


    long read_count = 0;
    unsigned long totalcap = 0;

    std::string bc1_word1;
    while (bc1_reader.getline(bc1_word1)) {

        std::vector<std::string> bc1_words;
        std::vector<std::string> bc2_words;
        std::vector<std::string> r1_words;
        std::vector<std::string> r2_words;

        bc1_words.push_back(std::move(bc1_word1));
        int lcount = 3;
        get_lines_to_vec(bc1_reader, bc1_words, lcount);
        lcount = 4;
        get_lines_to_vec(bc2_reader, bc2_words, lcount);
        get_lines_to_vec(reader1, r1_words, lcount);
        get_lines_to_vec(reader2, r2_words, lcount);

        ++read_count;

        // The entire short read consists bc1_index
        std::string bc1_local = bc1_words[1];
        std::string bc2_local = bc2_words[1];

        const auto& bc1_tuple = bc1_loader.match_barcode(bc1_local, mm_bc1);
        const auto& bc2_tuple = bc2_loader.match_barcode(bc2_local, mm_bc2);

        std::string barcode_str = write_entry(bc1_tuple, bc2_tuple);
        total_count++;
        //std::cout << "\n";

       // Get the stagger sequence using bc2

        unsigned int allcap = updateMaps(barcode_str, bc1_words, bc2_words, 
            r1_words, r2_words);
        totalcap += allcap;

        if (totalcap > total_allowed) {
            writeMapsToFile(count_log);
            // Write all the data in the respective files sequentially
            totalcap = 0;
        }

        
        if (read_count % 1000000 == 0) {
            std::cout << "read_count: " << read_count << "\n";
            std::cout << "count_map_exact_size: " << count_map_exact.size() << 
                 "\n";
            auto t = std::time(nullptr);
			auto tm = *std::localtime(&t);

		    std::ostringstream oss;
		    oss << std::put_time(&tm, "%d-%m-%Y %H-%M-%S");
		    auto str = oss.str();

		    std::cout << str << std::endl;
            //break;
        }
        if (read_count == total_run) {break;}
    }
    writeMapsToFile(count_log);
    std::string outfile_str = outdirpath + "/" + prefix_str + "_outfile.txt";
	
    write_outfile(outfile_str);
    std::cout << "Unique exact read count:\t" << unique_exact_count << "\n";
    std::cout << "Unique nonexact read count:\t" << unique_nonexact_count 
        << "\n";
    std::cout << "Ambiguous read count:\t" << ambiguous_count << "\n";
    std::cout << "No match read count:\t" << no_match_count << "\n";
    std::cout << "Total read count:\t" << total_count << "\n";
    unsigned long total_count_t = unique_exact_count + unique_nonexact_count + 
        ambiguous_count + no_match_count;
    assert(total_count_t == total_count);
    
}

bool dual_bc_splitter::has_suffix(const std::string &str, const std::string &suffix) {
       return str.size() >= suffix.size() &&
           str.compare(str.size() - suffix.size(), suffix.size(), suffix) == 0;
}


void dual_bc_splitter::compress_files() {

    std::string suffix(".fastq");
    for (auto& lfile : outfile_all_set) {
        if (has_suffix(lfile, suffix)) {
            // COnvert it to .gz file
            std::string cmd = "gzip " + lfile;
            std::cout << "Compressing: " << cmd << "\n";
            std::system(cmd.c_str());
        }
    }
}

int main(int argc, char* argv[]) { 

	dual_bc_splitter csm;

	
	bool all_set = true;
	try {
		all_set = csm.parse_args(argc, argv);	
	} catch(std::exception& e) {
		std::cerr << "error: " << e.what() << "\n";
        //lbs.print_help();
        return 1;

	} catch (...) {
		//lbs.print_help();
		return 0;
	}

 	if (!all_set) {
		csm.print_help();
		return 0;
	}

	csm.initialize();
	try {
        csm.load_barcodes();
		csm.core_engine();
        csm.compress_files();
	} catch(std::invalid_argument& e) {
        std::cerr << "error: " << e.what() << "\n";
		//lbs.print_help();
		return 1;
    }

	//csm.write_log();
        
    return 0;
}


