import janis as j
from janis.unix.tools.echo import Echo

w = j.WorkflowBuilder("test")

w.input("string1", j.String, default="Hello")
w.input("string2", j.String, default=" World")
w.input("not_echo_string2", j.Boolean, default=False)
w.step("echo_string1", Echo(inp=w.string1))
if not w.not_echo_string2.value:
    w.step("echo_string2", Echo(inp=w.string2))

if __name__ == "__main__":
    w.translate("cwl", to_disk=True)