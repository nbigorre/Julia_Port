module fortVar

export lkGet, fortSet, fortSetStr, fortGet, fortGetArr, fortSetArr, @lkGet, @fortSet, @fortGet, @fortSetArr, @fortGetArr

function fortSet(v_name, value, v_type)
    local ptr = cglobal(("header_mp_" * v_name * "_", "./PSOM_LIB.so"), v_type)
    unsafe_store!(ptr, value)
end

function fortSetStr(v_name, value, v_type, size)
    local ptr = cglobal(("header_mp_" * v_name * "_", "./PSOM_LIB.so"), v_type)
    #unsafe_store!(ptr, pointer(value))
    local arr = fill(UInt8(32), size)
    for i in eachindex(value)
        arr[i] = value[i]
    end
    Base.memcpy(ptr, pointer(arr), size)
end

function fortGet(v_name, v_type)
    local ptr = cglobal(("header_mp_" * v_name * "_", "./PSOM_LIB.so"), v_type)
    return unsafe_load(ptr)
end

function fortGetArr(v_name, v_type, dims)
    local ptr = cglobal(("header_mp_" * v_name * "_", "./PSOM_LIB.so"), v_type)
    return unsafe_wrap(Array{v_type}, ptr, dims)
end

function fortSetArr(v_name, value, v_type, size)
    local ptr = cglobal(("header_mp_" * v_name * "_", "./PSOM_LIB.so"), v_type)
    Base._memcpy!(ptr, value, size)
end

function lkGet(name, type)
    return cglobal(("header_mp_" * name * "_", "./PSOM_LIB.so"), type)
end

macro lkGet(name, type)
    local msg = "header_mp_" * name * "_"
    return :(cglobal(($msg, "./PSOM_LIB.so"), $(esc(type))))
end

macro fortSet(v_name, value, v_type)
    local msg = "header_mp_" * v_name * "_"
    return :(unsafe_store!(cglobal(($msg, "./PSOM_LIB.so"), $(esc(v_type))), $(esc(value))))
end

macro fortGet(v_name, v_type)
    local msg = "header_mp_" * v_name * "_"
    return :(unsafe_load(cglobal(($msg, "./PSOM_LIB.so"), $(esc(v_type)))))
end

macro fortGetArr(v_name, v_type, dims)
    local msg = "header_mp_" * v_name * "_"
    return :(unsafe_wrap(Array{$(esc(v_type))}, cglobal(($msg, "./PSOM_LIB.so"),$(esc(v_type))), $(esc(dims))))
end

macro fortSetArr(v_name, value, v_type, size)
    local msg = "header_mp_" * v_name * "_"
    return :(Base._memcpy!(cglobal(($msg, "./PSOM_LIB.so"),$(esc(v_type))), $(esc(value)), sizeof($(esc(v_type))) * $(esc(size))))
end


end



#=hostname = Vector{UInt8}(undef, 256)
const NI=24  
const NJ=24
hostname = "haha"
fortSetStr("dirout", hostname, Ptr{UInt8})
hehe = fortGet("dirout", Ptr{UInt8})
print(unsafe_string(hehe))
=#