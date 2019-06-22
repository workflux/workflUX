import pexpect
p = pexpect.spawn("bash")
p.sendline("FINISH_STRING=DONE")
p.sendline("cwltool --help")
p.sendline("EXITCODE=$(echo $?)")
# p.sendline("sleep 2")
p.sendline('echo "${FINISH_STRING} EXITCODE:${EXITCODE}:${FINISH_STRING}"')
p.expect("DONE:EXITCODE:.*:DONE")
p.before
p.after
int(p.after.decode().split(":")[2].strip())


