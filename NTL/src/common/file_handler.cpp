#include "file_handler.h"


namespace NTL::fh
{
    std::fstream json_to_file(const json& j, const std::string& name)
    {
        std::fstream stream;
        stream.open(name + ".json", std::ios_base::out);

        if (stream)
        {
            stream << j.dump();
        }
        else
            throw(std::runtime_error("Couldn't open " + name + " for writing"));

        stream.close();

        return stream;
    }

    json file_to_json(const std::string& name)
    {
        std::fstream stream;
        stream.open(name + ".json", std::ios_base::in);
        json j;

        if (stream)
        {
            j = json::parse(stream);
        }
        else
            throw(std::runtime_error("Couldn't open " + name + " for reading"));

        return j;
    }

    json ntl_to_json(const NTL& ntl, const std::string& readme)
    {
        if (ntl.get_Cn().empty())
            throw(std::logic_error("Attempted to write an empty NTL into a json object"));

        return {
            {"json_type", "NTL"},
            {"Z0", ntl.get_Z0()},
            {"e_r", ntl.get_er()},
            {"d", ntl.get_d()},
            {"Cn", ntl.get_Cn()},
            {"readme", readme}
        };
    }

    NTL json_to_ntl(const json& j)
    {
        if (j.at("json_type") != "NTL")
            throw(std::logic_error("Attempted to read an NTL from a different json object"));

        return NTL(j.at("Z0"), j.at("e_r"), j.at("d"), j.at("Cn"));
    }

    std::fstream ntl_to_file(const NTL& ntl, const std::string& name, const std::string& readme)
    {
        return json_to_file(ntl_to_json(ntl, readme), name);
    }

    NTL file_to_ntl(const std::string& name)
    {
        return json_to_ntl(file_to_json(name));
    }
}