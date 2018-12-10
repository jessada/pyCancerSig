#!/bin/bash

die () {
    echo >&2 "[exception] $@"
    echo >&2 "$usage"
    exit 1
}

msg_to_out ()
{
    local message="$1"
    echo "$message" 1>&2
    write_log "$message"
}

info_msg ()
{
    local message="${1:-}"

    INFO_MSG_FORMAT="## [INFO] `date "+%Y-%m-%d %H:%M:%S,%3N"` - %s"
    local formated_msg=`printf "$INFO_MSG_FORMAT" "$message"`
    msg_to_out "$formated_msg"
}

debug_msg ()
{
    local message="$1"

    DEBUG_MSG_FORMAT="## [DEBUG] `date "+%Y-%m-%d %H:%M:%S,%3N"` - %s"
    local formated_msg=`printf "$DEBUG_MSG_FORMAT" "$message"`
    if [ "$dev_mode" == "On" ]
    then
        msg_to_out "$formated_msg"
    else
        write_log "$formated_msg"
    fi
}

display_param ()
{
    local PARAM_PRINT_FORMAT="  %-40s%s"
    local param_name=$1
    local param_val=$2

    local msg=`printf "$PARAM_PRINT_FORMAT" "$param_name"":" "$param_val"`
    info_msg "$msg"
}

new_section_txt ()
{
    local section_message="$1"
    info_msg
    info_msg "************************************************** $section_message **************************************************"
}

new_sub_section_txt ()
{
    local sub_section_message="$1"
    info_msg
    info_msg ">>>>>>>>>>>>>>>>>>>> $sub_section_message <<<<<<<<<<<<<<<<<<<<"
}

eval_cmd ()
{
    cmd="$1"
    
    msg="executing: $cmd"

    info_msg
    info_msg "$msg"
    eval "$cmd"
}
