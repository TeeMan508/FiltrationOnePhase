#include "params.h"

#include <utility>
#include <sstream>
#include <fstream>
#include <string.h>


void Parameters::clear_comments(std::string &line) {
    if (line.find('#') == std::string::npos) return;
    line.erase(line.find('#'));
}

template<typename T, typename P>
T Parameters::remove_if(T beg, T end, P pred)
{
    T dest = beg;
    for (T itr = beg;itr != end; ++itr)
        if (!pred(*itr))
            *(dest++) = *itr;
    return dest;
}

std::pair<std::string, std::string> Parameters::parse_line(std::string line) {
    clear_comments(line);
    line.erase(remove_if(line.begin(), line.end(), isspace), line.end());

    if (line.find('=') !=  std::string::npos) {
        std::string name = line.substr(0, line.find('='));
        std::string value = line.substr(line.find('=') + 1);
        return std::pair(name, value);
    }

    if (line.find(':') !=  std::string::npos) {
        std::string name = line.substr(0, line.find(':'));
        return std::pair(name, "");
    }

    return std::pair("", "");
}

std::string Parameters::get_block_from_string(std::string sample) {
    std::stringstream value_stream(sample);
    std::string string_to;
    std::string accumulator;
    int counter = 0, line_idx = 0;

    while(std::getline(value_stream, string_to, '\n')) {
        std::string line = string_to;
        line.erase(remove_if(line.begin(), line.end(), isspace), line.end());
        if (line_idx > 0) accumulator += string_to + '\n';

        if (string_to.find(':') !=  std::string::npos) counter++;

        if (string_to.find('/') !=  std::string::npos
            and line.length() == 1) counter--;

        if (counter == 0) {
            return accumulator;
        }
        line_idx++;
    }
    return sample;
}

std::string Parameters::get_lines_from_number(std::string sample, int idx=0) {
    std::stringstream value_stream(sample);
    std::string string_to;
    int counter = 0;
    std::string accumulator;

    while(std::getline(value_stream, string_to, '\n')) {
        if (counter>=idx) accumulator += string_to + '\n';
        counter ++;

    }

    return accumulator;
}

int Parameters::measure_block_string(std::string block) {

    std::stringstream value_stream(block);
    std::string string_to;
    int counter = 0;
    std::string accumulator;
    while(std::getline(value_stream, string_to, '\n')) {
        counter ++;
    }
    return counter;
}

Parameters::Parameters(const std::string& sample) {
    form_from_string(*this, sample);
}

void Parameters::form_from_string(Parameters &params, std::string sample) {
    if (sample.empty()) return;

    std::stringstream value_stream(sample);
    std::string string_to;
    int current_line_idx = 0;
    int skip_lines = 0;
    while(std::getline(value_stream, string_to, '\n')) {
        std::pair<std::string, std::string> buffer = parse_line(string_to);
        if (skip_lines>0) {
            current_line_idx++;
            skip_lines--;
            continue;
        }
        if (!buffer.first.empty() and buffer.second.empty()) {
            if (params.blocks.find(buffer.first) != params.blocks.end()) {
                return;
            }
            std::string text_slice = get_lines_from_number(sample, current_line_idx);
            std::string block_string = get_block_from_string(text_slice);
            skip_lines = measure_block_string(block_string);
            params.blocks[buffer.first] = Parameters(block_string);


        }
        if (!buffer.first.empty() and !buffer.second.empty()) {
            params.fields[buffer.first] = buffer.second;
        }

        current_line_idx++;
    }
}

void Parameters::load(std::string file) {
    std::ifstream input_file_stream(file);
    std::stringstream buffer;
    std::string file_content, string_to;

    buffer << input_file_stream.rdbuf();
    while(std::getline(buffer, string_to, '\n')) {
        file_content += string_to + '\n';
    }
    form_from_string(*this, file_content);
}

void Parameters::print_tmp() {
    for (auto const&x : this->fields) {
        std::cout << x.first << " " << x.second << std::endl;
    }
    for (auto const&x : this->blocks) {
        std::cout << std::endl;
        std::cout << x.first << std::endl;
        for (auto const&y : x.second.fields) {
            std::cout << y.first << " " << y.second << std::endl;
        }
        for (auto const&z : x.second.blocks) {
            std::cout << std::endl;
            std::cout << z.first << std::endl;
            for (auto const&yz : z.second.fields) {
                std::cout << yz.first << " " << yz.second << std::endl;
            }
        }
    }
}

void Parameters::print_fields() {
    for (auto const&x : this->fields) {
        std::cout << "\t" << x.first << "=" << x.second << std::endl;
    }
}

int tab_counter = 1;
int recursion_steps = 0;
void Parameters::print_block(std::string block_name) {
    // tab_counter++;
    std::cout << std::endl;
    std::cout << std::string(tab_counter-1, '\t') << "Block: |" << block_name << '|' << std::endl;

    this->print_fields();

    for (auto &x : this->blocks) {
        x.second.print_block(x.first);
    }
}


void Parameters::print_all() {
    std::cout << "Fields:" << std::endl;
    this->print_fields();
    std::cout << "Blocks:" << std::endl;
    for (std::pair<Name, Parameters> x : this->blocks) {
        x.second.print_block(x.first);

    }
}

void Parameters::set(Name name, Value value) {

    if (this->fields.find(name) != this->fields.end()) {
        auto it = this->fields.find(name);
        it->second = value;
        return;
    }

    for (auto &x : this->blocks) {
        x.second.set(name, value);
    }

    this->fields[name] = value;
}


Value Parameters::get(Name name) {
    Value result;
    if (this->fields.find(name) != this->fields.end()) return this->fields[name];

    for (auto &x : this->blocks) {
        result = x.second.get(name);
    }

    return result;
}
//
// Value Parameters::get_or_set(Name name, Value value="") {
//     if (this->get(name).empty()) this->set(name, value);
//
//     return this->get(name);
// }

Value Parameters::operator[](Name const &name) {
    return this->get(name);
}



