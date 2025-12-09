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
            {"M", ntl.get_M()},
            {"readme", readme}
        };
    }

    NTL json_to_ntl(const json& j)
    {
        if (j.at("json_type") != "NTL")
            throw(std::logic_error("Attempted to read an NTL from a different json object"));

        return NTL(j.at("Z0"), j.at("e_r"), j.at("d"), j.at("Cn"), j.at("M"));
    }

    json wpd_to_json(const WPD::WPD& wpd, const std::string& readme)
    {
        if (wpd.get_ntl2().get_Cn().empty() || wpd.get_ntl3().get_Cn().empty() || !wpd.get_R())
            throw(std::logic_error("Attempted to write an empty NTL into a json object"));

        json j1 = {
            {"json_type", "WPD"},
            {"Z0", wpd.get_Z0()},
            {"e_r", wpd.get_er()},
            {"R", wpd.get_R()},
            {"readme", readme}
        };

        json j2 = ntl_to_json(wpd.get_ntl2()); 
        j2.erase("json_type");
        j2.erase("readme");
        j2.erase("Z0");
        j2.erase("e_r");        
        rename_json_element<double>(j2, "d", "d-2");
        rename_json_element<std::vector<double>>(j2, "Cn", "Cn2");
        rename_json_element<double>(j2, "M", "M-2");

        json j3 = ntl_to_json(wpd.get_ntl3()); 
        j3.erase("json_type");
        j3.erase("readme");
        j3.erase("Z0");
        j3.erase("e_r");
        rename_json_element<double>(j3, "d", "d-3");
        rename_json_element<std::vector<double>>(j3, "Cn", "Cn3");
        rename_json_element<double>(j3, "M", "M-3");

        json j_total = j1;         
        j_total.merge_patch(j2);   
        j_total.merge_patch(j3);

        return j_total;       
    }

    WPD::WPD json_to_wpd(const json& j)
    {
        if (j.at("json_type") != "WPD")
            throw(std::logic_error("Attempted to read a WPD from a different json object"));

        return WPD::WPD(NTL(j.at("Z0"), j.at("e_r"), j.at("d-2"), j.at("Cn2"), j.at("M-2")),
            NTL(j.at("Z0"), j.at("e_r"), j.at("d-3"), j.at("Cn3"), j.at("M-3")),
            j.at("R"));
    }

    std::fstream ntl_to_file(const NTL& ntl, const std::string& name, const std::string& readme)
    {
        return json_to_file(ntl_to_json(ntl, readme), name);
    }

    NTL file_to_ntl(const std::string& name)
    {
        return json_to_ntl(file_to_json(name));
    }

    std::fstream wpd_to_file(const WPD::WPD& wpd, const std::string& name, const std::string& readme)
    {
        return json_to_file(wpd_to_json(wpd, readme), name);
    }
    WPD::WPD file_to_wpd(const std::string& name)
    {
        return json_to_wpd(file_to_json(name));
    }
}