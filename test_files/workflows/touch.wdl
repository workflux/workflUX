version 1.0

workflow touch_wf {
  input {
    String filename
  }
  call touch {
      input:
        filename = filename
  }
}

task touch {
    input {
        String filename
    }
    command {
        touch "${filename}"
    }
    output {
        String out = read_string(stdout())
    }
}