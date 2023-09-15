# import subprocess

# # 执行 C++ 程序的命令
# cpp_command = "~/文档/uhd/host/build/examples/mytest --tx-freq  5000e6 --tx-rate 100e6 --rx-freq 5000e6 --rx-rate 100e6 --wave-type HELTF20 --wave-freq 2e6 --tx-gain 25 --rx-gain 25 --ampl 0.7 --settling 0.6 --type float --rx-channels 1 --tx-channels 1  --ref internal --tx-bw 20e6 --rx-bw 20e6 --tx-ant TX/RX --rx-ant RX2 --spb 28800  --file \"/mnt/tmp/usrp_20_multi.dat\""
# cpp_process = subprocess.Popen(cpp_command)

# # 执行 Python 程序的命令
# python_command = "/bin/python3 /home/dellx/文档/godirect-examples-main/python/gdx_getting_started_graphing.py"
# python_process = subprocess.Popen(python_command)

# # 等待子进程执行完毕
# cpp_process.wait()
# python_process.wait()
import subprocess
import threading

def run_cpp_program():
    cpp_command = "~/文档/uhd/host/build/examples/mytest \
    --tx-freq  5815e6 --tx-rate 200e6 --rx-freq 5815e6 --rx-rate 200e6 \
    --wave-type HELTF160 --wave-freq 2e6 --tx-gain 0 --rx-gain 30 --ampl 1.0 --settling 0.6 \
    --type float --rx-channels 1 --tx-channels 1  --ref internal --tx-bw 160e6 --rx-bw 160e6 \
    --tx-ant TX/RX --rx-ant RX2 --spb 230400 --file \"/mnt/tmp/usrp_160_0724_x2_multi5.dat\""
    subprocess.run(cpp_command, shell=True)

def run_python_program():
    python_command = "/bin/python3 \
    /home/dellx/文档/godirect-examples-main/python/gdx_getting_started_graphing.py \
    \"/mnt/tmp/usrp_160_0724_x2_multi5.csv\""
    subprocess.run(python_command, shell=True)

# 创建并启动线程
cpp_thread = threading.Thread(target=run_cpp_program)
cpp_thread.start()

python_thread = threading.Thread(target=run_python_program)
python_thread.start()

# 等待线程执行完毕
cpp_thread.join()
python_thread.join()

