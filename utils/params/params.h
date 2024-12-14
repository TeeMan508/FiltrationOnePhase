#ifndef PARAMS_H
#define PARAMS_H


#include <string>
#include <vector>
#include <iostream>
#include <map>

typedef std::string Name;
typedef std::string Value;
typedef std::map<Name, std::string> Fields;


class Parameters {
    Fields fields;
    std::map<Name, Parameters> blocks;

    // Utils functions
    static std::pair<Name, Value> parse_line(std::string line);
    static void clear_comments(std::string &line);

    template<typename T, typename P>
    static T remove_if(T beg, T end, P pred);

    static void form_from_string(Parameters &params, std::string sample);
    static std::string get_block_from_string(std::string sample);
    static std::string get_lines_from_number(std::string sample, int idx);
    static int measure_block_string(std::string block);

    Value search(Name);
public:
    // Constructors
    Parameters() = default;
    Parameters(const Parameters & b) = default;
    explicit Parameters(const std::string& value);

    // Methods
    void load(std::string filename);
    void print_tmp();
    void print_fields();
    void print_block(Name block_name);
    void print_all();
    Value get(Name name);
    void set(Name name, Value value);
    // Value get_or_set(Name name, Value value);

    //Operators
    Parameters& operator = (Parameters const& second) = default;


    Value operator [] (Name const& name);


};


#endif