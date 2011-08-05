/***********************************************************************
TCPSocket - Wrapper class for TCP sockets ensuring exception safety.

Copyright (c) 2003 The Regents of the University of California, through
Lawrence Berkeley National Laboratory, University of California at
Davis, and Lawrence Livermore National Laboratory, subject to any
required approvals from the U.S. Department of Energy.

This source code is part of the ProteinShop software.

ProteinShop is copyrighted and your use is under license, subject to
any required approvals from the U.S. Department of Energy.  For details
or questions, you may contact Berkeley Lab's Technology Transfer
Department at TTD@lbl.gov (Re:  ProteinShop; CR-1877)

NOTICE OF U.S. GOVERNMENT RIGHTS.  ProteinShop was developed under
funding from the U.S. Government which consequently retains certain
rights as follows: the U.S. Government has been granted for itself and
others acting on its behalf a paid-up, nonexclusive, irrevocable,
worldwide license in ProteinShop to reproduce, prepare derivative
works, and perform publicly and display publicly.  Beginning five (5)
years after the date permission to assert copyright is obtained from the
U.S. Department of Energy, and subject to any subsequent five (5) year
renewals, the U.S. Government is granted for itself and others acting on
its behalf a paid-up, nonexclusive, irrevocable, worldwide license in
ProteinShop to reproduce, prepare derivative works, distribute copies
to the public, perform publicly and display publicly, and to permit
others to do so.

Written by Oliver Kreylos.
***********************************************************************/

#ifndef TCPSOCKET_INCLUDED
#define TCPSOCKET_INCLUDED

#include <string>
#include <stdexcept>

class TCPSocket
{
    /* Embedded classes: */
    public:
    class PipeError:public std::runtime_error // Exception for unexpected connection termination
    {
        /* Constructors and destructors: */
        public:
        PipeError(const std::string& what_arg)
            :std::runtime_error(what_arg)
        {
        };
    };
    class TimeOut:public std::runtime_error // Exception for time-outs when waiting for data
    {
        /* Constructors and destructors: */
        public:
        TimeOut(const std::string& what_arg)
            :std::runtime_error(what_arg)
        {
        };
    };
    
    /* Elements: */
    private:
    int socketFd; // Internal socket file descriptor
    
    /* Constructors and destructors: */
    private:
    TCPSocket(int sSocketFd) // Creates a TCPSocket wrapper around an existing socket file descriptor (without copying)
        :socketFd(sSocketFd)
    {
    };
    public:
    TCPSocket(int portId,int backlog); // Creates a socket on the local host and starts listening; if portId is negative, random free port is assigned
    TCPSocket(const char* hostname,int portId); // Creates a socket connected to a remote host
    TCPSocket(const TCPSocket& source); // Copy constructor
    ~TCPSocket(void); // Closes a socket
    
    /* Methods: */
    int getFd(void) // Returns low-level socket file descriptor
    {
        return socketFd;
    };
    TCPSocket& operator=(const TCPSocket& source); // Assignment operator
    int getPortId(void) const; // Returns port ID assigned to a socket
    TCPSocket accept(void) const; // Waits for an incoming connection on a listening socket and returns a new socket connected to the initiator
    
    /* I/O methods: */
    bool waitForData(long timeoutSeconds,long timeoutMicroseconds,bool throwException =true) const; // Waits for incoming data on TCP socket; returns true if data is ready; (optionally) throws exception if wait times out
    size_t read(void* buffer,size_t count); // Reads raw buffer from TCP socket; returns number of bytes read
    void blockingRead(void* buffer,size_t count); // Reads raw buffer from TCP socket; blocks until data completely read
    void blockingWrite(const void* buffer,size_t count); // Writes raw buffer to TCP socket; blocks until data completely written
};

#endif
