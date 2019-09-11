from subprocess import Popen, PIPE
p = Popen('bash', stdin=PIPE, stdout=PIPE, bufsize=1)
cmd_string="TEST=HELLO"
p.stdin.write(bytes(cmd_string+'\n','utf-8'))
p.stdin.flush()
cmd_string="echo $TEST"
p.stdin.write(bytes(cmd_string+'\n','utf-8'))
p.stdin.flush()
p.stdout.readline()
p.terminate()
p.poll()
