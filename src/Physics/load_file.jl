using HDF5

function iterate_named_tuple(ps, path, file, keys_file, keys_ps)
    for (name, value) in pairs(ps)
        new_path = vcat(path, name)
        if value isa NamedTuple
            keys_ps = iterate_named_tuple(value, new_path, file, keys_file, keys_ps)
        else
            if keys_ps % 2 == 0
                value .= file[keys_file[keys_ps-1]]
            else
                tmp = Array(file[keys_file[keys_ps+1]])
                value .= tmp'
                
            end
            keys_ps += 1
        end
    end

    return keys_ps
end

function OpenHDF5(filename, ps)
    h5open(filename, "r") do file
        keys_file = keys(file)
        keys_ps = 1 
        iterate_named_tuple(ps, Symbol[], file, keys_file, keys_ps)  
    end
end
