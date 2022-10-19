# Extract and output scan information from minicbf files

# Assumptions:
# 1. files are named in some form of xxxx_scanno(_)frameno.cbf
# 2. frames are sequential
# 3. The first frame will be scan 1 frame 1
# 4. filenames are the same length
# 5. axis names are drawn from the lists below
# 6. "<axis>_increment" signals the increment

using ArgParse

# Configuration information

const rot_axes = ("chi", "phi", "detector_2theta", "two_theta", "omega")
const trans_axes = ("detector_distance", "dx", "trans")

is_trans_axis(axis) = axis in trans_axes
is_rot_axis(axis) = axis in rot_axes

get_scan_info(frame_dir) = begin

    frame_catcher, all_names = get_scan_frame_fmt(frame_dir)

    # index them all for speed

    all_frames = map(all_names) do an
        m = match(frame_catcher, an)
        (m["scan"],parse(Int, m["frame"])) => an
    end

    all_frames = Dict(all_frames)
    
    all_scans = unique!(map( x -> x[1], collect( keys( all_frames ) ) ) )
    println("$(length(all_scans)) scans found")
    scan_info = Dict()

    axes = union(rot_axes, trans_axes)
    
    for s in all_scans
        
        snames = filter( x-> x[1] == s, keys(all_frames))
        frames = [ x[2] for x in snames ]

        @debug "$(length(frames)) frames found for scan $s"

        # Get information for first and last frame

        fname = joinpath( frame_dir, all_frames[(s, 1)])
        vals, scan_ax, scan_incr, et, wl = get_frame_info(fname, axes)
        start = vals[scan_ax]

        fname = joinpath( frame_dir, all_frames[(s, length(frames))])
        vals, _, _, _, _ = get_frame_info(fname, axes)
        finish = vals[scan_ax]

        # Check increment and range match
        
        if !isapprox(start + scan_incr * (length(frames)-1), finish, atol=1e-6)
            throw(error("Scan range does not match increment: $start to $finish, $(length(frames)-1) steps of $scan_incr"))
        end
       
        details = Dict("frames" => length(frames),
                       "axis" => scan_ax,
                       "incr" => scan_incr,
                       "time" => et,
                       "start" => start,
                       "range" => scan_incr * length(frames),
                       "wavelength" => wl)
        
        scan_info[s] = (vals, details)
    end

    return scan_info, all_frames
    
end

# Deduce scan/frame naming convention
get_scan_frame_fmt(frame_dir) = begin
    
    all_names = readdir(frame_dir)
    filter!(x-> x[end-3:end] == ".cbf", all_names)
    sort!(all_names)

    rr = r".*_(?<scan>[0-9]+1)(?<sep>0|_)(?<frame>[0-9]+1)\.cbf"
    m = match(rr,all_names[1])

    s_pos = m.offsets[1]
    s_len = length(m["scan"])
    f_len = length(m["frame"])
    if m["sep"] == "0"
        f_pos = m.offsets[2]
        f_len = f_len + 1
    else
        f_pos = m.offsets[3]
    end

    fr_sc_regex = ".*_(?<scan>[0-9]{$s_len})_?(?<frame>[0-9]{$f_len})"

    @assert !isnothing(match(Regex(fr_sc_regex), all_names[end]))
    
    return Regex(fr_sc_regex), all_names
end

"""
Return any values found for provided axes. All axes converted to lowercase.
"""
get_frame_info(fname, axes) = begin
    
    header = readuntil(fname, "--CIF-BINARY-FORMAT-SECTION--")
    lines = lowercase.(split(header, "\n"))

    ax_vals = map(ax -> (lowercase(ax), (get_header_value(lines, ax))), axes)
    filter!( x-> x[2] != nothing, ax_vals)
    ax_vals = convert_units.(ax_vals)

    ax_incr = map(ax -> (lowercase(ax), get_header_value(lines, ax*"_increment")), axes)
    filter!( x-> x[2] != nothing, ax_incr)
    ax_incr = convert_units.(ax_incr)
    
    et, _ = get_header_value(lines, "exposure_time")
    wl, _ = get_header_value(lines, "wavelength")

    scan_ax = findfirst( x -> !isapprox(x[2], 0, atol = 1e-6), ax_incr)
    
    return Dict(ax_vals), ax_incr[scan_ax][1], ax_incr[scan_ax][2], et, wl
    
end

"""
    Get the value following the string given in matcher and units if present
"""
get_header_value(lines, matcher) = begin

    rr = Regex("$matcher[ =]+")
    one_line = filter( x-> !isnothing(match(rr, x)), lines)
    if length(one_line) != 1
        return nothing
    end

    one_line = one_line[]
    #@debug "Extracting from" one_line

    m = match(Regex("$matcher[ =]+(?<val>[A-Za-z0-9+-.]+) +(?<units>[A-Za-z.]+)"), one_line)
    val = strip(m["val"])
    units = strip(m["units"])

    #@debug "To get value" val

    return parse(Float64, val), units
end

"""
   Detect any non-mm translations and convert
"""
convert_units(ax_val) = begin
        name, (val, units) = ax_val
        if units == "m"
            val = val * 1000
        elseif units == "cm"
            val = val * 10
        end
        name, val
end

"""
Axes in `scan_info` and `axis_dict` should already be lowercase.
"""
rename_axes!(scan_info, axis_dict) = begin

    if length(axis_dict) == 0 return end

    for (s, v) in scan_info

        vals, dets = v
        for (a, x) in vals

            if a in keys(axis_dict)
                vals[axis_dict[a]] = x
                delete!(vals, a)
            end
        end

        if dets["axis"] in keys(axis_dict)
            dets["axis"] = axis_dict[dets["axis"]]
        end
    end
    
end

#============= Output routines ========================#

"""
We identify each frame by its sequence number overall (not just within
its own scan.)
"""

output_scan_info(scan_info, all_frames, output_file, new_url; prepend_dir = "", arch = nothing) = begin

    sl = create_scan_list(scan_info)
    exp_info = Dict([s=>scan_info[s][2]["time"] for s in keys(scan_info)])

    op = isnothing( output_file ) ? stdout : open(output_file, "w")

    generate_wavelength(op, scan_info)
    generate_scan_settings(op, scan_info)
    generate_scan_info(op, sl)
    generate_step_info(op, sl, exp_info)
    generate_array_info(op, sl)
    generate_ids(op, sl)
    generate_external_ids(op, new_url, all_frames, sl, prepend_dir, comp=arch)
    
end

create_scan_list(scan_info) = begin

    # Create scan list of (scanid, frame_no) where
    # frame_no is just the number of frames

    scans = sort(collect(keys(scan_info)))
    slist = [(s, scan_info[s][2]["frames"]) for s in scans]

    @debug "Scan list" slist
    
    return slist
end

generate_wavelength(op, scan_info; rad_type = "xray") = begin
    wl = first(scan_info).second[2]["wavelength"]
    println(op,"_diffrn_radiation_wavelength.id 1")
    println(op,"_diffrn_radiation_wavelength.value $wl")
    println(op,"_diffrn_radiation.type $rad_type")
end

"""

Information we have obtained from the cbf file

        details = Dict("frames"
                       "axis"
                       "incr"
                       "time"
                       "start"
                       "range")
        

"""
generate_scan_settings(op, scan_info) = begin
    header = """loop_
_diffrn_scan_axis.scan_id
_diffrn_scan_axis.axis_id
_diffrn_scan_axis.displacement_start
_diffrn_scan_axis.displacement_increment
_diffrn_scan_axis.displacement_range
_diffrn_scan_axis.angle_start
_diffrn_scan_axis.angle_increment
_diffrn_scan_axis.angle_range
"""
    println(op, header)
    for s in sort(collect(keys(scan_info)))
        
        axes, dets = scan_info[s]
        
        for (a, val) in axes

            step, range = 0, 0
            if a == dets["axis"]
                step = dets["incr"]
                range = dets["range"]
            end
            
            if is_trans_axis(a)
                println(op, "SCAN$s  $a $val $step $range . . .")
            else
                println(op, "SCAN$s  $a . . . $val $step $range")
            end
        end
        println(op, "")

    end

    
end

"""

Array_id is just IMAGE (a single detector module). The
binary_id is incremented with frame, and all of them
are located externally so external_id is also incremented
together with binary_id.
"""
generate_ids(op, scan_list) = begin
    header = """loop_
_array_data.array_id
_array_data.binary_id
_array_data.external_data_id
"""
    println(op, header)
    ctr = 0
    for (s,f) in scan_list
        for i in 1:f
            ctr += 1
            println(op, "   IMAGE $ctr $ctr")
        end
    end
    println(op,"")
end

"""

`fulluri` is the location of a single archive file. The individual files
on the local storage are assumed to be at the same locations relative to
this top-level directory.
"""
generate_external_ids(op, fulluri, all_frames, scan_list, prepend_dir; comp="TBZ",fmt="CBF") = begin
    header = """loop_
_array_data_external_data.id
_array_data_external_data.format
_array_data_external_data.uri
"""
    if !isnothing(arch)
        header = header *
"""_array_data_external_data.archive_format
_array_data_external_data.archive_path"""
    end

    # If comp is nothing, then each frame has a separate URI and is treated as a local file.
    # If comp is something, then each frame is relative to a single URI and is output
    # separately.
    println(op, header)
    ctr = 0
    for (s,f) in scan_list
        for i in 1:f
            ctr += 1
            fname = all_frames[(s,i)]
            print(op, "  $ctr $fmt $fulluri")

            # A too-clever-by-half way of optionally live constructing a URL
            
            if !isnothing(arch)
                print(op, "  $comp $prepend_dir")
            end
            
            println(op, "/" * "$fname")
        end
    end
    println(op,"")
end

"""
Fill in the scan information. We number the frames from
the start
"""
generate_scan_info(op, scan_list) = begin
    header = """loop_
_diffrn_scan.id 
_diffrn_scan.frame_id_start
_diffrn_scan.frame_id_end
_diffrn_scan.frames"""
    println(op, header)
    stptr = 1
    for (s,f) in scan_list
        endptr = stptr + f - 1
        println(op, "SCAN$s    frm$stptr  frm$endptr  $f")
        stptr = endptr + 1
    end
    println(op,"")
end

"""
Fill in information about steps
"""
generate_step_info(op, scan_list, scan_times) = begin
    header = """loop_
_diffrn_scan_frame.frame_id
_diffrn_scan_frame.scan_id
_diffrn_scan_frame.frame_number
_diffrn_scan_frame.integration_time"""
    println(op, header)
    ctr = 0
    for (s,n) in scan_list
        for f in 1:n
            ctr += 1
            println(op, "frm$ctr   SCAN$s    $f $(scan_times[s])")
        end
    end
    println(op, "")
end

generate_array_info(op, scan_list) = begin
    header = """loop_
_diffrn_data_frame.id
_diffrn_data_frame.detector_element_id
_diffrn_data_frame.array_id
_diffrn_data_frame.binary_id"""
    ctr = 0
    println(op, header)
    for (s,n) in scan_list
        for f in 1:n
            ctr += 1
            println(op, "frm$ctr  ELEMENT  IMAGE $(ctr)")
        end
    end
end

determine_archive(location) = begin

    arch = nothing
    
    if location == []
        location = "file://" * frame_dir
        
    else
        location = parsed_args["location"][]
        long_end = location[end-6:end]
        short_end = location[end-2:end]
        
        if short_end == "tgz" || long_end == ".tar.gz"
            arch = "TGZ"
        elseif short_end == "tbz" || long_end == "tar.bz2"
            arch = "TBZ"
        elseif short_end == "zip"
            arch = "ZIP"
        end

    end

    return arch, location
        
end

parse_cmdline(d) = begin

    s = ArgParseSettings(d)

    @add_arg_table! s begin
        "-l", "--location"
        help = "Final URL of files in output."
        nargs = 1
        "-i", "--include"
        help = "Include directory name as part of frame location info in output"
        nargs = 0
        "-o", "--output"
        help = "Output file to write to, if missing stdout"
        nargs = 1
        "-a", "--axis"
        nargs = 2
        metavar = ["cbf", "new"]
        action = "append_arg"
        help = "Change axis name from <cbf> to <new> in output file to match goniometer axis definitions. May be used multiple times for multiple axis renaming" 
        "directory"
        help = "Directory containing scan frames in minicbf format"
        required = true
    end

    parse_args(s)

end

if abspath(PROGRAM_FILE) == @__FILE__

    # Process arguments
    
    parsed_args = parse_cmdline("Extract scan information from minicbf files")
    frame_dir = parsed_args["directory"]
    do_output = parsed_args["output"] != nothing
    prepend_dir = parsed_args["include"] ? splitpath(frame_dir)[end] : ""

    # Analyse CBF files
    
    scan_info, all_frames = get_scan_info(frame_dir)

    # Rename axes

    @debug "axis" parsed_args["axis"]
    if length(parsed_args["axis"]) > 0
        lower_axes = Dict([lowercase(x[1]) => x[2] for x in parsed_args["axis"]])
        @debug "lower" lower_axes
        rename_axes!(scan_info, lower_axes)
    end
    
    # Output CIF fragment

    arch, location = determine_archive(parsed_args["location"])
    
    out_file = parsed_args["output"] != [] ? parsed_args["output"][] : nothing
    if do_output
        output_scan_info(scan_info, all_frames, out_file, location,
                         prepend_dir = prepend_dir,
                         arch = arch)
    end
end
