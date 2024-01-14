#include <iostream>
#include <cstring>
#include <sys/socket.h>
#include <netinet/in.h>
#include <arpa/inet.h>
#include <unistd.h>

int main() {
    int serverSocket, clientSocket;
    struct sockaddr_in serverAddr, clientAddr;
    socklen_t addrSize;

    // Create a socket
    serverSocket = socket(AF_INET, SOCK_STREAM, 0);
    if (serverSocket < 0) {
        std::cerr << "Error in socket creation\n";
        return -1;
    }

    // Set up the server address struct
    memset(&serverAddr, 0, sizeof(serverAddr));
    serverAddr.sin_family = AF_INET;
    serverAddr.sin_addr.s_addr = INADDR_ANY;
    serverAddr.sin_port = htons(8080); // Port number

    // Bind the socket to the IP and port
    if (bind(serverSocket, (struct sockaddr *)&serverAddr, sizeof(serverAddr)) < 0) {
        std::cerr << "Error in binding\n";
        return -1;
    }

    // Listen for incoming connections
    listen(serverSocket, 5);

    // Accept a client connection
    addrSize = sizeof(clientAddr);
    clientSocket = accept(serverSocket, (struct sockaddr *)&clientAddr, &addrSize);
    if (clientSocket < 0) {
        std::cerr << "Error in accepting connection\n";
        return -1;
    }

    float buffer[4096];
    size_t totalBytesReceived = 0;

    while (totalBytesReceived < sizeof(buffer)) {
        ssize_t bytesRead = recv(clientSocket, reinterpret_cast<char*>(buffer) + totalBytesReceived, sizeof(buffer) - totalBytesReceived, 0);
        
        if (bytesRead <= 0) {
            // 接收出错或连接关闭
            break;
        }

        totalBytesReceived += bytesRead;
    }

    size_t numFloatsReceived = totalBytesReceived / sizeof(float);
    // 处理接收到的浮点数数据
    for (size_t i = 0; i < numFloatsReceived; ++i) {
        std::cout << "Received float[" << i << "]: " << buffer[i] << std::endl;
    }

    // Close the sockets
    close(clientSocket);
    close(serverSocket);

    return 0;
}
