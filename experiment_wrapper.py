import sys, os, argparse, subprocess, time
from concurrent.futures import ThreadPoolExecutor
from threading import Condition


def prepare_args(argv):
    for i in range(len(argv)):
        if argv[i] == '--':
            command_with_args = argv[i+1:]
            args_to_tag_output = argv[1:i]
            return args_to_tag_output, command_with_args
    return None, None


def make_argparser():
    parser = argparse.ArgumentParser(description="Tag standard output/error")
    parser.add_argument("-o", "--output-file", default="", help="Specify file to print output to (default: stdout)")
    parser.add_argument("-e", "--error-file", default="", help="Specify file to print error to (default: where stdout ends up)")
    parser.add_argument("--tee", dest='tee', default=False, action='store_true', help="If output is redirected to a file, still print it to stdout as well")
    parser.add_argument("--time-format", default="[% 10.4f] ", help="Specify an alternative format string to tag lines with")
    parser.add_argument("--time-limit", type=float, default=-1, help="Specify a time limit for the command, enforced by killing the process (default: no timelimit)")
    parser.add_argument("--print-end-message", dest='print_end_message', default=False, action='store_true', help="Print a message to stdout as if the process printed it when the process ended")
    parser.add_argument("--dont-tag", dest="dont_tag", default=False, action='store_true', help="Do not tag or change the output at all.")
    return parser


def use_session_and_pgkill():
    return os.name != 'nt'


try:
    # import package that we only need on some platforms
    import signal
except ImportError:
    pass


def handle_commandline_args():
    args_to_us, command = prepare_args(sys.argv)
    if args_to_us is None or not command:
        print(f"Usage: {sys.argv[0]} [ARGUMENTS TO OUTPUT TAGGER] -- COMMAND [COMMAND ARGUMENTS]")
        make_argparser().print_help()
        sys.exit(1)
    return make_argparser().parse_args(args_to_us), command


class ProcessHandler:
    def __init__(self, args, command):
        self.args = args
        self.tag_format = self.args.time_format.encode('utf-8')
        self.command = command
        self.begin_time = None
        self.done_condition = Condition()
        self.command_is_done = False
        self.time_limit_expired = False
        self.process = None
        self.stdout_outputs = []
        self.stderr_outputs = []

    def __passed_time(self):
        return (time.time_ns() - self.begin_time) / 1000000000.0

    def __process_done(self):
        with self.done_condition:
            self.command_is_done = True
            self.done_condition.notify_all()

    def __timelimit_expired(self):
        self.time_limit_expired = True
        if self.process is not None:
            if not use_session_and_pgkill():
                self.process.kill()
            else:
                os.killpg(os.getpgid(self.process.pid), signal.SIGKILL)

    def __await_time_limit(self):
        with self.done_condition:
            while not self.command_is_done:
                passed = self.__passed_time()
                remaining_limit = 1000000.0 if self.args.time_limit < 0.0 else self.args.time_limit - passed
                if remaining_limit <= 0:
                    self.__timelimit_expired()
                    return
                self.done_condition.wait(remaining_limit)

    def __create_outputs(self):
        if self.args.output_file:
            self.stdout_outputs.append(open(self.args.output_file, "wb"))
            if self.args.tee:
                self.stdout_outputs.append(sys.stdout.buffer)
        else:
            self.stdout_outputs.append(sys.stdout.buffer)

        if self.args.error_file:
            ef = self.stdout_outputs if self.args.output_file == self.args.error_file else open(self.args.error_file, "wb")
            self.stderr_outputs.append(ef)
            if self.args.tee:
                self.stderr_outputs.append(sys.stderr.buffer)
        else:
            if self.args.output_file:
                self.stderr_outputs = self.stdout_outputs
            else:
                self.stderr_outputs = [sys.stderr.buffer]

    def __get_outputs(self, is_err):
        if is_err:
            return self.stderr_outputs
        else:
            return self.stdout_outputs

    def __process_stream_transparent(self, stream, is_err):
        outputs = self.__get_outputs(is_err)
        while True:
            data = stream.read1(16384)
            if not data:
                return
            for output in outputs:
                output.write(data)
                if output is sys.stdout.buffer or output is sys.stderr.buffer:
                    output.flush()

    def __write_tag(self, outputs):
        t = time.time_ns()
        tnow = (t - self.begin_time) / 1000000000.0
        tag = self.tag_format % tnow
        for output in outputs:
            output.write(tag)

    def __process_stream_tagging(self, stream, is_err):
        outputs = self.__get_outputs(is_err)
        last_was_newline = True
        while True:
            data = stream.read1(16384)
            if not data:
                return
            if last_was_newline:
                self.__write_tag(outputs)
                last_was_newline = (data[-1] == b'\n'[0])
            lines = data.splitlines()
            first = True
            for line in lines:
                if not first:
                    self.__write_tag(outputs)
                first = False
                for output in outputs:
                    output.write(line)
                    output.write(b'\n')
            for output in outputs:
                if output is sys.stdout.buffer or output is sys.stderr.buffer:
                    output.flush()
                
    def __process_stream(self, stream, is_err):
        if self.args.dont_tag:
            self.__process_stream_transparent(stream, is_err)
        else:
            self.__process_stream_tagging(stream, is_err)

    def __print_end_message(self):
        if self.time_limit_expired:
            for output in self.stdout_outputs:
                output.write(b'\n')
        self.__write_tag(self.stdout_outputs)
        for output in self.stdout_outputs:
            if self.time_limit_expired:
                output.write(b"PROCESS KILLED BY TIMEOUT\n")
            else:
                output.write(b"PROCESS ENDED\n")

    def run_process(self):
        self.begin_time = time.time_ns()
        self.process = subprocess.Popen(self.command, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=False,
                                        start_new_session=use_session_and_pgkill())
        try:
            self.__create_outputs()
            with ThreadPoolExecutor(4) as executor:
                r1 = executor.submit(self.__process_stream, self.process.stdout, False)
                r2 = executor.submit(self.__process_stream, self.process.stderr, True)
                r3 = executor.submit(self.__await_time_limit)
                self.process.wait()
                r1.result()
                r2.result()
                self.__process_done()
                r3.result()
            if self.args.print_end_message:
                self.__print_end_message()
            return self.process.returncode
        finally:
            for x in self.stdout_outputs:
                if x is sys.stdout.buffer or x is sys.stderr.buffer: continue
                x.close()
            for x in self.stderr_outputs:
                if x is sys.stdout.buffer or x is sys.stderr.buffer: continue
                if self.stdout_outputs and x is self.stdout_outputs[0]: continue
                x.close()


if __name__ == "__main__":
    args, command = handle_commandline_args()
    proc = ProcessHandler(args, command)
    res = proc.run_process()
    sys.exit(res)

