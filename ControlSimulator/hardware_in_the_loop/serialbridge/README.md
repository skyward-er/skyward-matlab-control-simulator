This README refers to the **serialbridge** module used for the Hardware In the Loop (HIL) simulation.

**Files**:
- `serialib.h` and `serialib.cpp`: used in order to build the C++ code for matlab
- `serialbridge.cpp`: the file that implements the serialbridge feature in C++
- `build.m`: used for building the serialbridge.cpp file
- `serialbridge.mexw64`: the compiled mex file that provides the serialbridge
    function for matlab (this is the only file needed for the simulation)

**Usage** of serialbridge in matlab:
- `serialbridge("Open", string_serialPort, uint_baudrate)`:
  - **"Open"**: specifies we want to open the serial port
  - string_serialPort: is the serial port we want to open (eg: "COM6")
  - uint_baudrate: is the baudrate of the port (eg: 256000)
- `serialbridge("Write", singleArray_Data)`:
  - **"Write"**: specifies that we want to write to the serial port
  - singleArray_Data: is the array of singles we want to write on serial (eg: [1 2 3 4.5 5.4])
- `singleArray_Data = serialbridge("Read", uint_nData)`;
  - **"Read"**: specifies that we want to read from the serial port
  - uint_nData: How many floats to read from serial (eg: 1)
  - singleArray_Data: array of floats read from the serial (eg: actuatorData)
- `serialbridge("Close", string_serialPort, uint_baudrate)`:
  - **"Close"**: specifies we want to close the serial port

**Example** in Matlab:
```
serialbridge("Open", "COM6", 256000); % Opens the serial port
serialbridge("Write", [1 2 3 4]);     % Sends the array "[1 2 3 4]" to the serial device
data = serialbridge("Read", 2);       % Receives 2 floats and stores them in the variable "data"
serialbridge("Close");                % Closes the serial port
```

Procedure in order to use the **MatlabTransceiver** module: https://git.skywarder.eu/scs/hermes/r2a-obsw/-/blob/Serial4Simulations-dev/src/tests/hardware_in_the_loop/README.md

***WARNING***:
* It's possible to open just ONE serial port with this module
* The files serialib.h and serialib.cpp were modified in order to fix compiling errors on windows
