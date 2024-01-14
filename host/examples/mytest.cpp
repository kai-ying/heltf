//
// Copyright 2010-2012,2014-2015 Ettus Research LLC
// Copyright 2018 Ettus Research, a National Instruments Company
//
// SPDX-License-Identifier: GPL-3.0-or-later
//
#include "heltf_CP_wavetable.hpp"
#include "heltf_wavetable.hpp"
#include <uhd/exception.hpp>
#include <uhd/types/tune_request.hpp>
#include <uhd/usrp/multi_usrp.hpp>
#include <uhd/utils/safe_main.hpp>
#include <uhd/utils/static.hpp>
#include <uhd/utils/thread.hpp>
#include <boost/algorithm/string.hpp>
#include <boost/bind.hpp>
#include <boost/filesystem.hpp>
#include <boost/format.hpp>
#include <boost/math/special_functions/round.hpp>
#include <boost/program_options.hpp>
#include <boost/thread/thread.hpp>
#include <csignal>
#include <fstream>
#include <iostream>
#include <cstring>
#include <sys/socket.h>
#include <netinet/in.h>
#include <arpa/inet.h>
#include <unistd.h>
#include <atomic>



namespace po = boost::program_options;
// std::mutex mtx;
// size_t index = 0;
boost::mutex indexMutex; // 声明一个互斥量

boost::condition_variable tranmitCP_condition;
boost::mutex cv_m;
std::atomic<bool> restart_transmit(false);

/***********************************************************************
 * Signal handlers
 **********************************************************************/
static bool stop_signal_called = false;
void sig_int_handler(int)
{
    stop_signal_called = true;
}

/***********************************************************************
 * Utilities
 **********************************************************************/
//! Change to filename, e.g. from usrp_samples.dat to usrp_samples.00.dat,
//  but only if multiple names are to be generated.
std::string generate_out_filename(
    const std::string& base_fn, size_t n_names, size_t this_name)
{
    if (n_names == 1) {
        return base_fn;
    }

    boost::filesystem::path base_fn_fp(base_fn);
    base_fn_fp.replace_extension(boost::filesystem::path(
        str(boost::format("%02d%s") % this_name % base_fn_fp.extension().string())));
    return base_fn_fp.string();
}

<<<<<<< HEAD
/***********************************************************************
 * server_worker function
 * A function to be used as a boost::thread_group thread for Server Socket 
 **********************************************************************/
void server_work(size_t &index) {
    // 创建 Socket
    int serverSocket = socket(AF_INET, SOCK_STREAM, 0);
    if (serverSocket == -1) {
        std::cerr << "Error: Could not create socket\n";
        // return 1;
    }

    // 绑定地址和端口
    struct sockaddr_in serverAddr;
    serverAddr.sin_family = AF_INET;
    serverAddr.sin_port = htons(8080); // 端口号
    serverAddr.sin_addr.s_addr = INADDR_ANY;

    if (bind(serverSocket, (struct sockaddr *)&serverAddr, sizeof(serverAddr)) == -1) {
        std::cerr << "Error: Binding failed\n";
        close(serverSocket);
        // return 1;
    }

    // 监听连接
    if (listen(serverSocket, 5) == -1) {
        std::cerr << "Error: Listening failed\n";
        close(serverSocket);
        // return 1;
    }

    std::cout << "Server waiting for incoming connection...\n";

    // 接受客户端连接
    struct sockaddr_in clientAddr;
    socklen_t clientAddrSize = sizeof(clientAddr);
    int clientSocket = accept(serverSocket, (struct sockaddr *)&clientAddr, &clientAddrSize);
    if (clientSocket == -1) {
        std::cerr << "Error: Accepting client failed\n";
        close(serverSocket);
        // return 1;
    }

    std::cout << "Client connected\n";

    // 循环接收客户端发送的整数数据
    while (not stop_signal_called) {
        int receivedData;
        if (recv(clientSocket, &receivedData, sizeof(receivedData), 0) == -1) {
            std::cerr << "Error: Receiving data failed\n";
            close(clientSocket);
            close(serverSocket);
            // return 1;
        }

        std::cout << "Received data from client: " << receivedData << std::endl;
        if (abs(receivedData + 256) < 50 && receivedData != 0) {
            // boost::lock_guard<boost::mutex> lock(indexMutex); // 在访问 index 之前加锁
            // ADD SOCKET Message to Client
            const char* message = "OK!";
            send(clientSocket, message, strlen(message), 0);
            std::cout << "Synchronization has done!" << std::endl;
            //////////////////////////////////////////////////////////////
            // LKY DO IT: receive change CP message
            char CP_flag[2]; // 设置接收缓冲区
            memset(CP_flag, 0, sizeof(CP_flag)); // 初始化缓冲区
            if (recv(clientSocket, CP_flag, sizeof(CP_flag), 0) == -1) {
                std::cerr << "Error: Receiving data failed\n";
                close(clientSocket);
                // return 1;
            }
            if (strcmp(CP_flag, "CP") == 0) {
                // TO DO IT: change to CP data
                std::cout << "Need to change to data with CP Now!" << std::endl;
                {
                    // 上锁
                    boost::lock_guard<boost::mutex> lock(cv_m);
                    // 重置重启标志位
                    restart_transmit.store(true);
                }
                // 条件变量提醒线程
                tranmitCP_condition.notify_one();
            }

            boost::this_thread::sleep(boost::posix_time::seconds(100));
        }
        else {
            // // 通过锁定互斥量来修改共享数据
            // std::cout << "Initial index: " << index << std::endl;
            boost::lock_guard<boost::mutex> lock(indexMutex); // 在访问 index 之前加锁
            index -= receivedData+256;
            // std::cout << "Index After compensation: " << index << std::endl;

            // Sleeping for 0.1 second
            boost::this_thread::sleep(boost::posix_time::milliseconds(100));
        }
    }

    // 关闭连接
    close(clientSocket);
    close(serverSocket);
}

/***********************************************************************
 * transmitCP_worker function
 * A function to be used as a boost::thread_group thread for transmitting
 **********************************************************************/
void transmitCP_worker(std::vector<std::complex<float>> CPbuff,
    CPwave_table_class CPwave_table,
    uhd::tx_streamer::sptr tx_streamer,
    uhd::tx_metadata_t metadata,
    size_t step,
    size_t &index,
    size_t cnt,
    int num_channels)
{
    std::vector<std::complex<float>*> CPbuffs(num_channels, &CPbuff.front());
    {
        boost::unique_lock<boost::mutex> lock(cv_m);
        // 等待条件变量通知
        tranmitCP_condition.wait(lock);
        std::cout << "transmitCP_worker start!!" << std:: endl;
    }
    boost::this_thread::sleep(boost::posix_time::microseconds(150));
    // send data until the signal handler gets called
    while (not stop_signal_called && restart_transmit) {
        if (++cnt == 12) {
            cnt = 0;
            index--;
        }
        // fill the buffer with the waveform
        for (size_t n = 0; n < CPbuff.size(); n++) {
            // buff[n] = wave_table();
            // if (not restart_transmit) {
            //     if (index % 2048 == 0) index += 256;
            // }
            CPbuff[n] = CPwave_table(index += step);
            // std::cout << buff[n] << std::endl;
            // std::cout << index << std:: endl;
        }
        // std::cout << buff.size() << std::endl;

        // send the entire contents of the buffer
        tx_streamer->send(CPbuffs, CPbuff.size(), metadata);

        metadata.start_of_burst = false;
        metadata.has_time_spec  = false;
    }
    std::cout << std::endl << "transmitCP_worker finished!!!" << std::endl << std::endl;

    // send a mini EOB packet
    metadata.end_of_burst = true;
    tx_streamer->send("", 0, metadata);
}

=======
std::vector<std::complex<float>> tx_bb;
>>>>>>> db0a79b7223892f6974b3c4ce4f6c9373df2501b

/***********************************************************************
 * transmit_worker function
 * A function to be used as a boost::thread_group thread for transmitting
 **********************************************************************/
void transmit_worker(std::vector<std::complex<float>> buff,
    wave_table_class wave_table,
    uhd::tx_streamer::sptr tx_streamer,
    uhd::tx_metadata_t metadata,
    size_t step,
    size_t &index,
    size_t cnt,
    int num_channels)
{
    std::vector<std::complex<float>*> buffs(num_channels, &buff.front());
    
    // if (not restart_transmit) {
    //     // std::cout << std::endl << "Now transmit_worker!!!" << std::endl << std::endl;
    
    // } else { // restart transmmit_thread
    //     transmit_thread.interrupt_all(); //中断之前的发送线程
    //     transmit_thread.join_all();  // 等待线程结束
    //     boost::thread_group transmit_thread;
    //     transmit_thread.create_thread(boost::bind(
    //         &transmitCP_worker, buff, CPwave_table, tx_stream, md, step, boost::ref(index), cnt, num_channels));
    //     restart_transmit = false;
    //     std::cout << std::endl << "Have Changed to transmitCP_worker!!!" << std::endl << std::endl;
    // }


    // send data until the signal handler gets called
    while (not stop_signal_called && not restart_transmit) {
        // Tx_wavefrom index
        // std::cout << "TX_wavefrom index" << index <<std::endl;
        // boost::lock_guard<boost::mutex> lock(indexMutex); // 在访问 index 之前加锁
        if (++cnt == 12) {
            cnt = 0;
            index--;
        }
        // fill the buffer with the waveform
        for (size_t n = 0; n < buff.size(); n++) {
<<<<<<< HEAD
            // buff[n] = wave_table();
            // if (not restart_transmit) {
            //     if (index % 2048 == 0) index += 256;
            // }
            buff[n] = wave_table(index += step);
=======

            // index = (index + step) % wave_table_len;
            // buff[n] = wave_table(index);
            buff[n] = wave_table();
            // if (tx_bb.size() < 2304) tx_bb.push_back(buff[n]);
>>>>>>> db0a79b7223892f6974b3c4ce4f6c9373df2501b
            // std::cout << buff[n] << std::endl;
            // std::cout << index << std:: endl;
        }        

        // send the entire contents of the buffer
        tx_streamer->send(buffs, buff.size(), metadata);

        metadata.start_of_burst = false;
        metadata.has_time_spec  = false;
    }

    // // //send bb signal
    // for (size_t n = 0; n < 2304; n++) {
    //     // std::cout << buff.size() << std::endl;
    //     std::cout << buff[n] << std::endl;
    // }
    // for (size_t j = 0; j < 2304; j++) {
    //     std::cout << tx_bb[j] << std::endl;
    // }


    // send a mini EOB packet
    // metadata.end_of_burst = true;
    // tx_streamer->send("", 0, metadata);
    std::cout << std::endl << "transmit_worker finished!!!" << std::endl << std::endl;
    // boost::this_thread::sleep(boost::posix_time::seconds(200));
    // // transmit_thread.interrupt_all(); //中断之前的发送线程
    // transmit_thread.join_all();  // 等待线程结束
    // boost::thread_group transmit_thread;
    // transmit_thread.create_thread(boost::bind(
    //     &transmitCP_worker, buff, CPwave_table, tx_stream, md, step, boost::ref(index), cnt, num_channels));
    // restart_transmit = false;
    // std::cout << std::endl << "Have Changed to transmitCP_worker!!!" << std::endl << std::endl;
    
}


/***********************************************************************
 * recv_to_file function
 **********************************************************************/
template <typename samp_type>
void recv_to_file(uhd::usrp::multi_usrp::sptr usrp,
    const std::string& cpu_format,
    const std::string& wire_format,
    const std::string& file,
    size_t samps_per_buff,
    int num_requested_samples, // int / long?
    double settling_time,
    std::vector<size_t> rx_channel_nums)
{

    std::vector<samp_type> rx_data;

    int num_total_samps = 0;
    // create a receive streamer
    uhd::stream_args_t stream_args(cpu_format, wire_format);
    stream_args.channels             = rx_channel_nums;
    uhd::rx_streamer::sptr rx_stream = usrp->get_rx_stream(stream_args);

    // Prepare buffers for received samples and metadata
    uhd::rx_metadata_t md;
    std::vector<std::vector<samp_type>> buffs(
        rx_channel_nums.size(), std::vector<samp_type>(samps_per_buff));
    // create a vector of pointers to point to each of the channel buffers
    std::vector<samp_type*> buff_ptrs;
    for (size_t i = 0; i < buffs.size(); i++) {
        buff_ptrs.push_back(&buffs[i].front());

    }

    // Create one ofstream object per channel
    // (use shared_ptr because ofstream is non-copyable)
    std::vector<boost::shared_ptr<std::ofstream>> outfiles;
    for (size_t i = 0; i < buffs.size(); i++) {
        const std::string this_filename = generate_out_filename(file, buffs.size(), i);
        outfiles.push_back(boost::shared_ptr<std::ofstream>(
            new std::ofstream(this_filename.c_str(), std::ofstream::binary)));
    }


    UHD_ASSERT_THROW(outfiles.size() == buffs.size());
    UHD_ASSERT_THROW(buffs.size() == rx_channel_nums.size());
    bool overflow_message = true;
    double timeout =
        settling_time + 0.1f; // expected settling time + padding for first recv

    // setup streaming
    uhd::stream_cmd_t stream_cmd((num_requested_samples == 0)
                                     ? uhd::stream_cmd_t::STREAM_MODE_START_CONTINUOUS
                                     : uhd::stream_cmd_t::STREAM_MODE_NUM_SAMPS_AND_DONE);
    stream_cmd.num_samps  = num_requested_samples;
    stream_cmd.stream_now = false;
    stream_cmd.time_spec  = uhd::time_spec_t(settling_time);
    rx_stream->issue_stream_cmd(stream_cmd);

    while (not stop_signal_called
           and (num_requested_samples > num_total_samps or num_requested_samples == 0)) {
        size_t num_rx_samps = rx_stream->recv(buff_ptrs, samps_per_buff, md, timeout);
        timeout             = 0.1f; // small timeout for subsequent recv

        if (md.error_code == uhd::rx_metadata_t::ERROR_CODE_TIMEOUT) {
            std::cout << boost::format("Timeout while streaming") << std::endl;
            break;
        }
        if (md.error_code == uhd::rx_metadata_t::ERROR_CODE_OVERFLOW) {
            if (overflow_message) {
                overflow_message = false;
                std::cerr
                    << boost::format(
                           "Got an overflow indication. Please consider the following:\n"
                           "  Your write medium must sustain a rate of %fMB/s.\n"
                           "  Dropped samples will not be written to the file.\n"
                           "  Please modify this examp0722pear again.\n")
                           % (usrp->get_rx_rate() * sizeof(samp_type) / 1e6);
            }
            continue;
        }
        if (md.error_code != uhd::rx_metadata_t::ERROR_CODE_NONE) {
            throw std::runtime_error(
                str(boost::format("Receiver error %s") % md.strerror()));
        }
<<<<<<< HEAD

        num_total_samps += num_rx_samps;

        // for (size_t i = 0; i < outfiles.size(); i++) {
        //     outfiles[i]->write(
        //         (const char*)buff_ptrs[i], num_rx_samps / 100 * sizeof(samp_type));
        // }
=======
        num_total_samps += num_rx_samps;
        // 不降采样 /50
        for (size_t i = 0; i < outfiles.size(); i++) {
            outfiles[i]->write(
                (const char*)buff_ptrs[i], num_rx_samps  * sizeof(samp_type) / 100); // / 100 
        }

>>>>>>> db0a79b7223892f6974b3c4ce4f6c9373df2501b
    }
    
    //recv complex float signal
    for (size_t i = 0; i < outfiles.size(); i++) {
        // std::cout << outfiles.size() << std::endl;
        // std::cout << sizeof(samp_type) << std::endl; 

        for (int t = 0; t < 2304; t++) {
            rx_data.push_back(buff_ptrs[i][t]);
            std::cout << rx_data[t] << std::endl;
            // std::cout << rx_data.size() << std::endl;
            // std::cout << buff_ptrs[i][t] << std::endl;
        }
    }


    // Shut down receiver
    stream_cmd.stream_mode = uhd::stream_cmd_t::STREAM_MODE_STOP_CONTINUOUS;
    rx_stream->issue_stream_cmd(stream_cmd);


    // Close files
    for (size_t i = 0; i < outfiles.size(); i++) {
        outfiles[i]->close();
    }

}


/***********************************************************************
 * Main function
 **********************************************************************/
int UHD_SAFE_MAIN(int argc, char* argv[])
{
    // transmit variables to be set by po
    std::string tx_args, wave_type, tx_ant, tx_subdev, ref, otw, tx_channels;
    double tx_rate, tx_freq, tx_gain, wave_freq, tx_bw;
    float ampl;

    // receive variables to be set by po
    std::string rx_args, file, type, rx_ant, rx_subdev, rx_channels;
    size_t total_num_samps, spb, CPspb;
    double rx_rate, rx_freq, rx_gain, rx_bw;
    double settling;

    // setup the program options
    po::options_description desc("Allowed options");
    // clang-format off
    desc.add_options()
        ("help", "help message")
        ("tx-args", po::value<std::string>(&tx_args)->default_value(""), "uhd transmit device address args")
        ("rx-args", po::value<std::string>(&rx_args)->default_value(""), "uhd receive device address args")
        ("file", po::value<std::string>(&file)->default_value("usrp_samples.dat"), "name of the file to write binary samples to")
        ("type", po::value<std::string>(&type)->default_value("short"), "sample type in file: double, float, or short")
        ("nsamps", po::value<size_t>(&total_num_samps)->default_value(0), "total number of samples to receive")
        ("settling", po::value<double>(&settling)->default_value(double(0.2)), "settling time (seconds) before receiving")
        ("spb", po::value<size_t>(&spb)->default_value(0), "samples per buffer, 0 for default")
        ("CPspb", po::value<size_t>(&CPspb)->default_value(0), "samples per buffer, 0 for default")

        ("tx-rate", po::value<double>(&tx_rate), "rate of transmit outgoing samples")
        ("rx-rate", po::value<double>(&rx_rate), "rate of receive incoming samples")
        ("tx-freq", po::value<double>(&tx_freq), "transmit RF center frequency in Hz")
        ("rx-freq", po::value<double>(&rx_freq), "receive RF center frequency in Hz")
        ("ampl", po::value<float>(&ampl)->default_value(float(0.3)), "amplitude of the waveform [0 to 0.7]")
        ("tx-gain", po::value<double>(&tx_gain), "gain for the transmit RF chain")
        ("rx-gain", po::value<double>(&rx_gain), "gain for the receive RF chain")
        ("tx-ant", po::value<std::string>(&tx_ant), "transmit antenna selection")
        ("rx-ant", po::value<std::string>(&rx_ant), "receive antenna selection")
        ("tx-subdev", po::value<std::string>(&tx_subdev), "transmit subdevice specification")
        ("rx-subdev", po::value<std::string>(&rx_subdev), "receive subdevice specification")
        ("tx-bw", po::value<double>(&tx_bw), "analog transmit filter bandwidth in Hz")
        ("rx-bw", po::value<double>(&rx_bw), "analog receive filter bandwidth in Hz")
        // ("wave-type", po::value<std::string>(&wave_type)->default_value("SINE"), "waveform type (HELTF20, HELTF40, HELTF80, HELTF160)")
        ("wave-type", po::value<std::string>(&wave_type)->default_value("HELTF20"), "waveform type (HELTF20, HELTF40, HELTF80, HELTF160)")
        ("wave-freq", po::value<double>(&wave_freq)->default_value(0), "waveform frequency in Hz")
        ("ref", po::value<std::string>(&ref)->default_value("internal"), "clock reference (internal, external, mimo)")
        ("otw", po::value<std::string>(&otw)->default_value("sc16"), "specify the over-the-wire sample mode")
        ("tx-channels", po::value<std::string>(&tx_channels)->default_value("0"), "which TX channel(s) to use (specify \"0\", \"1\", \"0,1\", etc)")
        ("rx-channels", po::value<std::string>(&rx_channels)->default_value("0"), "which RX channel(s) to use (specify \"0\", \"1\", \"0,1\", etc)")
        ("tx-int-n", "tune USRP TX with integer-N tuning")
        ("rx-int-n", "tune USRP RX with integer-N tuning")
    ;
    // clang-format on
    po::variables_map vm;
    po::store(po::parse_command_line(argc, argv, desc), vm);
    po::notify(vm);

    // print the help message
    if (vm.count("help")) {
        std::cout << boost::format("UHD TXRX Loopback to File %s") % desc << std::endl;
        return ~0;
    }

    // create a usrp device
    std::cout << std::endl;
    std::cout << boost::format("Creating the transmit usrp device with: %s...") % tx_args
              << std::endl;
    uhd::usrp::multi_usrp::sptr tx_usrp = uhd::usrp::multi_usrp::make(tx_args);
    std::cout << std::endl;
    std::cout << boost::format("Creating the receive usrp device with: %s...") % rx_args
              << std::endl;
    uhd::usrp::multi_usrp::sptr rx_usrp = uhd::usrp::multi_usrp::make(rx_args);

    // always select the subdevice first, the channel mapping affects the other settings
    if (vm.count("tx-subdev"))
        tx_usrp->set_tx_subdev_spec(tx_subdev);
    if (vm.count("rx-subdev"))
        rx_usrp->set_rx_subdev_spec(rx_subdev);

    // detect which channels to use
    std::vector<std::string> tx_channel_strings;
    std::vector<size_t> tx_channel_nums;
    boost::split(tx_channel_strings, tx_channels, boost::is_any_of("\"',"));
    for (size_t ch = 0; ch < tx_channel_strings.size(); ch++) {
        size_t chan = std::stoi(tx_channel_strings[ch]);
        if (chan >= tx_usrp->get_tx_num_channels()) {
            throw std::runtime_error("Invalid TX channel(s) specified.");
        } else
            tx_channel_nums.push_back(std::stoi(tx_channel_strings[ch]));
    }
    std::vector<std::string> rx_channel_strings;
    std::vector<size_t> rx_channel_nums;
    boost::split(rx_channel_strings, rx_channels, boost::is_any_of("\"',"));
    for (size_t ch = 0; ch < rx_channel_strings.size(); ch++) {
        size_t chan = std::stoi(rx_channel_strings[ch]);
        if (chan >= rx_usrp->get_rx_num_channels()) {
            throw std::runtime_error("Invalid RX channel(s) specified.");
        } else
            rx_channel_nums.push_back(std::stoi(rx_channel_strings[ch]));
    }

    // Lock mboard clocks
    if (vm.count("ref")) {
        tx_usrp->set_clock_source(ref);
        rx_usrp->set_clock_source(ref);
    }

    std::cout << boost::format("Using TX Device: %s") % tx_usrp->get_pp_string()
              << std::endl;
    std::cout << boost::format("Using RX Device: %s") % rx_usrp->get_pp_string()
              << std::endl;

    // set the transmit sample rate
    if (not vm.count("tx-rate")) {
        std::cerr << "Please specify the transmit sample rate with --tx-rate"
                  << std::endl;
        return ~0;
    }
    std::cout << boost::format("Setting TX Rate: %f Msps...") % (tx_rate / 1e6)
              << std::endl;
    tx_usrp->set_tx_rate(tx_rate);
    std::cout << boost::format("Actual TX Rate: %f Msps...")
                     % (tx_usrp->get_tx_rate() / 1e6)
              << std::endl
              << std::endl;

    // set the receive sample rate
    if (not vm.count("rx-rate")) {
        std::cerr << "Please specify the sample rate with --rx-rate" << std::endl;
        return ~0;
    }
    std::cout << boost::format("Setting RX Rate: %f Msps...") % (rx_rate / 1e6)
              << std::endl;
    rx_usrp->set_rx_rate(rx_rate);
    std::cout << boost::format("Actual RX Rate: %f Msps...")
                     % (rx_usrp->get_rx_rate() / 1e6)
              << std::endl
              << std::endl;

    // set the transmit center frequency
    if (not vm.count("tx-freq")) {
        std::cerr << "Please specify the transmit center frequency with --tx-freq"
                  << std::endl;
        return ~0;
    }

    for (size_t ch = 0; ch < tx_channel_nums.size(); ch++) {
        size_t channel = tx_channel_nums[ch];
        if (tx_channel_nums.size() > 1) {
            std::cout << "Configuring TX Channel " << channel << std::endl;
        }
        std::cout << boost::format("Setting TX Freq: %f MHz...") % (tx_freq / 1e6)
                  << std::endl;
        uhd::tune_request_t tx_tune_request(tx_freq);
        if (vm.count("tx-int-n"))
            tx_tune_request.args = uhd::device_addr_t("mode_n=integer");
        tx_usrp->set_tx_freq(tx_tune_request, channel);
        std::cout << boost::format("Actual TX Freq: %f MHz...")
                         % (tx_usrp->get_tx_freq(channel) / 1e6)
                  << std::endl
                  << std::endl;

        // set the rf gain
        if (vm.count("tx-gain")) {
            std::cout << boost::format("Setting TX Gain: %f dB...") % tx_gain
                      << std::endl;
            tx_usrp->set_tx_gain(tx_gain, channel);
            std::cout << boost::format("Actual TX Gain: %f dB...")
                             % tx_usrp->get_tx_gain(channel)
                      << std::endl
                      << std::endl;
        }

        // set the analog frontend filter bandwidth
        if (vm.count("tx-bw")) {
            std::cout << boost::format("Setting TX Bandwidth: %f MHz...") % (tx_bw / 1e6)
                      << std::endl;
            tx_usrp->set_tx_bandwidth(tx_bw, channel);
            std::cout << boost::format("Actual TX Bandwidth: %f MHz...")
                             % (tx_usrp->get_tx_bandwidth(channel) / 1e6)
                      << std::endl
                      << std::endl;
        }

        // set the antenna
        if (vm.count("tx-ant"))
            tx_usrp->set_tx_antenna(tx_ant, channel);
    }

    for (size_t ch = 0; ch < rx_channel_nums.size(); ch++) {
        size_t channel = rx_channel_nums[ch];
        if (rx_channel_nums.size() > 1) {
            std::cout << "Configuring RX Channel " << channel << std::endl;
        }

        // set the receive center frequency
        if (not vm.count("rx-freq")) {
            std::cerr << "Please specify the center frequency with --rx-freq"
                      << std::endl;
            return ~0;
        }
        std::cout << boost::format("Setting RX Freq: %f MHz...") % (rx_freq / 1e6)
                  << std::endl;
        uhd::tune_request_t rx_tune_request(rx_freq);
        if (vm.count("rx-int-n"))
            rx_tune_request.args = uhd::device_addr_t("mode_n=integer");
        rx_usrp->set_rx_freq(rx_tune_request, channel);
        std::cout << boost::format("Actual RX Freq: %f MHz...")
                         % (rx_usrp->get_rx_freq(channel) / 1e6)
                  << std::endl
                  << std::endl;

        // set the receive rf gain
        if (vm.count("rx-gain")) {
            std::cout << boost::format("Setting RX Gain: %f dB...") % rx_gain
                      << std::endl;
            rx_usrp->set_rx_gain(rx_gain, channel);
            std::cout << boost::format("Actual RX Gain: %f dB...")
                             % rx_usrp->get_rx_gain(channel)
                      << std::endl
                      << std::endl;
        }

        // set the receive analog frontend filter bandwidth
        if (vm.count("rx-bw")) {
            std::cout << boost::format("Setting RX Bandwidth: %f MHz...") % (rx_bw / 1e6)
                      << std::endl;
            rx_usrp->set_rx_bandwidth(rx_bw, channel);
            std::cout << boost::format("Actual RX Bandwidth: %f MHz...")
                             % (rx_usrp->get_rx_bandwidth(channel) / 1e6)
                      << std::endl
                      << std::endl;
        }

        // set the receive antenna
        if (vm.count("rx-ant"))
            rx_usrp->set_rx_antenna(rx_ant, channel);
    }

    // // for the const wave, set the wave freq for small samples per period
    // if (wave_freq == 0 and wave_type == "CONST") {
    //     wave_freq = tx_usrp->get_tx_rate() / 2;
    // }

    // error when the waveform is not possible to generate
    if (std::abs(wave_freq) > tx_usrp->get_tx_rate() / 2) {
        throw std::runtime_error("wave freq out of Nyquist zone");
    }
    if (tx_usrp->get_tx_rate() / std::abs(wave_freq) > wave_table_len / 2) {
        throw std::runtime_error("wave freq too small for table");
    }

    // pre-compute the waveform values
    const wave_table_class wave_table(wave_type, ampl);
    const size_t step = 1;
        // boost::math::iround(wave_freq / tx_usrp->get_tx_rate() * wave_table_len);
   
    // pre-compute the CPwaveform values
    const CPwave_table_class CPwave_table(wave_type, ampl);

    size_t index = 0;
    size_t cnt = 0;
    // index = 0;

    // create a transmit streamer
    // linearly map channels (index0 = channel0, index1 = channel1, ...)
    uhd::stream_args_t stream_args("fc32", otw);
    stream_args.channels             = tx_channel_nums;
    uhd::tx_streamer::sptr tx_stream = tx_usrp->get_tx_stream(stream_args);

    // allocate a buffer which we re-use for each channel
    if (spb == 0){
        spb = tx_stream->get_max_num_samps() * 10;
    }
    if (CPspb == 0){
        CPspb = tx_stream->get_max_num_samps() * 10;
    }
    std::cout << "spb: " << spb << std::endl;
    std::cout << "CPspb: " << CPspb << std::endl;
    std::cout << tx_stream->get_max_num_samps() << std::endl;  

    std::vector<std::complex<float>> buff(spb);
    int num_channels = tx_channel_nums.size();

    // LKY DO IT!
    std::vector<std::complex<float>> CPbuff(CPspb);


    // setup the metadata flags
    uhd::tx_metadata_t md;
    md.start_of_burst = true;
    md.end_of_burst   = false;
    md.has_time_spec  = true;
    md.time_spec = uhd::time_spec_t(0.5); // give us 0.5 seconds to fill the tx buffers

    // Check Ref and LO Lock detect
    std::vector<std::string> tx_sensor_names, rx_sensor_names;
    tx_sensor_names = tx_usrp->get_tx_sensor_names(0);
    if (std::find(tx_sensor_names.begin(), tx_sensor_names.end(), "lo_locked")
        != tx_sensor_names.end()) {
        uhd::sensor_value_t lo_locked = tx_usrp->get_tx_sensor("lo_locked", 0);
        std::cout << boost::format("Checking TX: %s ...") % lo_locked.to_pp_string()
                  << std::endl;
        UHD_ASSERT_THROW(lo_locked.to_bool());
    }
    rx_sensor_names = rx_usrp->get_rx_sensor_names(0);
    if (std::find(rx_sensor_names.begin(), rx_sensor_names.end(), "lo_locked")
        != rx_sensor_names.end()) {
        uhd::sensor_value_t lo_locked = rx_usrp->get_rx_sensor("lo_locked", 0);
        std::cout << boost::format("Checking RX: %s ...") % lo_locked.to_pp_string()
                  << std::endl;
        UHD_ASSERT_THROW(lo_locked.to_bool());
    }

    tx_sensor_names = tx_usrp->get_mboard_sensor_names(0);
    if ((ref == "mimo")
        and (std::find(tx_sensor_names.begin(), tx_sensor_names.end(), "mimo_locked")
                != tx_sensor_names.end())) {
        uhd::sensor_value_t mimo_locked = tx_usrp->get_mboard_sensor("mimo_locked", 0);
        std::cout << boost::format("Checking TX: %s ...") % mimo_locked.to_pp_string()
                  << std::endl;
        UHD_ASSERT_THROW(mimo_locked.to_bool());
    }
    if ((ref == "external")
        and (std::find(tx_sensor_names.begin(), tx_sensor_names.end(), "ref_locked")
                != tx_sensor_names.end())) {
        uhd::sensor_value_t ref_locked = tx_usrp->get_mboard_sensor("ref_locked", 0);
        std::cout << boost::format("Checking TX: %s ...") % ref_locked.to_pp_string()
                  << std::endl;
        UHD_ASSERT_THROW(ref_locked.to_bool());
    }

    rx_sensor_names = rx_usrp->get_mboard_sensor_names(0);
    if ((ref == "mimo")
        and (std::find(rx_sensor_names.begin(), rx_sensor_names.end(), "mimo_locked")
                != rx_sensor_names.end())) {
        uhd::sensor_value_t mimo_locked = rx_usrp->get_mboard_sensor("mimo_locked", 0);
        std::cout << boost::format("Checking RX: %s ...") % mimo_locked.to_pp_string()
                  << std::endl;
        UHD_ASSERT_THROW(mimo_locked.to_bool());
    }
    if ((ref == "external")
        and (std::find(rx_sensor_names.begin(), rx_sensor_names.end(), "ref_locked")
                != rx_sensor_names.end())) {
        uhd::sensor_value_t ref_locked = rx_usrp->get_mboard_sensor("ref_locked", 0);
        std::cout << boost::format("Checking RX: %s ...") % ref_locked.to_pp_string()
                  << std::endl;
        UHD_ASSERT_THROW(ref_locked.to_bool());
    }

    if (total_num_samps == 0) {
        std::signal(SIGINT, &sig_int_handler);
        std::cout << "Press Ctrl + C to stop streaming..." << std::endl;
    }

    // reset usrp time to prepare for transmit/receive
    std::cout << boost::format("Setting device timestamp to 0...") << std::endl;
    tx_usrp->set_time_now(uhd::time_spec_t(0.0));

    // start transmit worker thread
    boost::thread_group transmit_thread;
    transmit_thread.create_thread(boost::bind(
        &transmit_worker, buff, wave_table, tx_stream, md, step, boost::ref(index), cnt, num_channels));
    
    // LKY DO！
    transmit_thread.create_thread(boost::bind(
        &transmitCP_worker, CPbuff, CPwave_table, tx_stream, md, step, boost::ref(index), cnt, num_channels));

    // Create a boost::thread_group
    boost::thread_group server_Group;
    // Start the thread
    server_Group.create_thread(boost::bind(&server_work, boost::ref(index)));


    // // recv to file
    // if (type == "double")
    //     recv_to_file<std::complex<double>>(
    //         rx_usrp, "fc64", otw, file, spb, total_num_samps, settling, rx_channel_nums);
    // else if (type == "float")
    //     recv_to_file<std::complex<float>>(
    //         rx_usrp, "fc32", otw, file, spb, total_num_samps, settling, rx_channel_nums);
    // else if (type == "short")
    //     recv_to_file<std::complex<short>>(
    //         rx_usrp, "sc16", otw, file, spb, total_num_samps, settling, rx_channel_nums);
    // else {
    //     // clean up transmit worker
    //     stop_signal_called = true;
    //     transmit_thread.join_all();
    //     throw std::runtime_error("Unknown type " + type);
    // }
    if (type == "float") {
        while (not stop_signal_called) {

        };
    }

    // clean up transmit worker
    stop_signal_called = true;
    transmit_thread.join_all();

    // finished
    std::cout << std::endl << "Done!" << std::endl << std::endl;
    return EXIT_SUCCESS;
}
