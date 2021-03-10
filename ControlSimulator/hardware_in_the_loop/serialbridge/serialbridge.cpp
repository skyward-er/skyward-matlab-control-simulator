/* Copyright (c) 2020 Skyward Experimental Rocketry
 * Authors: Luca Erbetta
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in
 * all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.  IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
 * THE SOFTWARE.
 */

#include <string>

#include "mex.hpp"
#include "mexAdapter.hpp"

/* [WARNING]
 * I had to modify the serialib.h header and the serialib.cpp cpp file in order
 * to get rid of the compile errors on windows */
#include "serialib.h"

using namespace std;

static constexpr int INDEX_FUN         = 0;
static constexpr int INDEX_OPEN_PORT   = 1;
static constexpr int INDEX_OPEN_BAUD   = 2;
static constexpr int INDEX_READ_PARAM  = 1;
static constexpr int INDEX_WRITE_PARAM = 1;

static serialib serial;
static bool opened;

/**
 * @brief Matlab function to read/write from serial
 *
 * 1st parameter takes the action to execute:
 * - serialbridge("Open", string_serialPort, uint_baudrate):
 *   - "Open": specifies we want to open the serial port
 *   - string_serialPort: is the serial port we want to open (eg: "COM6")
 *   - uint_baudrate: is the baudrate of the port (eg: 256000)
 * - serialbridge("Write", singleArray_Data):
 *   - "Write": specifies that we want to write to the serial port
 *   - singleArray_Data: is the array of singles we want to write on serial
 * (eg:[1 2.2 3.33])
 * - singleArray_Data = serialbridge("Read", uint_nData);
 *   - "Read": specifies that we want to read from the serial port
 *   - uint_nData: How many floats to read from serial (eg: 1)
 *   - singleArray_Data: array of floats read from the serial (eg: actuatorData)
 * - serialbridge("Close", string_serialPort, uint_baudrate):
 *   - "Close": specifies we want to close the serial port
 *
 * eg, in Matlab:
 * serialbridge("Open", "COM6", 256000); % Opens the serial port
 *
 * serialbridge("Write", [1 2 3 4]);     % Sends the array "[1 2 3 4]" to the
 * serial device
 *
 * data = serialbridge("Read", 2);       % Receives 2 floats and stores them in
 * the variable "data"
 *
 * serialbridge("Close");                % Closes the serial port
 */
class MexFunction : public matlab::mex::Function
{
public:
    /**
     * @brief Entrypoint of the serialbridge matlab function
     */
    void operator()(matlab::mex::ArgumentList outputs,
                    matlab::mex::ArgumentList inputs)
    {
        std::shared_ptr<matlab::engine::MATLABEngine> matlabPtr = getEngine();
        matlab::data::ArrayFactory factory;

        checkArguments(outputs, inputs);

        if (inputs[INDEX_FUN][0] == "Open")
        {
            openSerial(inputs);
        }
        else if (inputs[INDEX_FUN][0] == "Close")
        {
            closeSerial(inputs);
        }
        else if (inputs[INDEX_FUN][0] == "Read")
        {
            read(outputs, inputs);
        }
        else if (inputs[INDEX_FUN][0] == "Write")
        {
            write(inputs);
        }
    }

    /**
     * @brief function used in order to open the serial.
     *
     * On matlab called with serialbridge("Open", portname, baudrate)
     * @param portname the port we want to open for serial communication
     * @param baudrate the baudrate of the port
     */
    void openSerial(matlab::mex::ArgumentList inputs)
    {
        std::shared_ptr<matlab::engine::MATLABEngine> matlabPtr = getEngine();
        matlab::data::ArrayFactory factory;

        string port = inputs[INDEX_OPEN_PORT][0];
        int baud    = inputs[INDEX_OPEN_BAUD][0];

        if (opened)
        {
            serial.closeDevice();
        }

        if (serial.openDevice(port.c_str(), baud) != 1)
        {
            opened = true;
            matlabPtr->feval(
                u"error", 0,
                std::vector<matlab::data::Array>(
                    {factory.createScalar("Error opening serial port.")}));
        }
        else
        {
            serial.flushReceiver();
            opened = true;

            disp("Serial port open.");
        }
    }

    /**
     * @brief function used in order to close the serial.
     *
     * On matlab called with serialbridge("Close")
     */
    void closeSerial(matlab::mex::ArgumentList inputs)
    {
        if (opened)
        {
            serial.closeDevice();
            opened = false;
            disp("Serial port closed.");
        }
        else
        {
            error("Serial port already closed.");
        }
    }

    /**
     * @brief function used in order to write on the serial.
     *
     * On matlab called with serialbridge("Write", data)
     * @param data array of singles to write on serial
     */
    void write(matlab::mex::ArgumentList &inputs)
    {
        if (!opened)
        {
            error("Serial port is not open!");
        }

        size_t numout = inputs[INDEX_WRITE_PARAM].getNumberOfElements();

        float *dataout = new float[numout];
        for (int i = 0; i < numout; i++)
        {
            dataout[i] = inputs[INDEX_WRITE_PARAM][i];
        }

        serial.writeBytes(dataout, sizeof(float) * numout);
    }

    /**
     * @brief function used in order to read on the serial.
     *
     * On matlab called with data = serialbridge("Read", n_floats)
     * @param n_floats int representing how many floats to read
     * @param data array of singles containing the data read from serial
     */
    void read(matlab::mex::ArgumentList &outputs,
              matlab::mex::ArgumentList &inputs)
    {
        if (!opened)
        {
            error("Serial port is not open!");
        }

        size_t numin    = inputs[INDEX_READ_PARAM][0];
        uint8_t *datain = new uint8_t[sizeof(float) * numin];
        size_t numbytes = 0;

        do
        {
            numbytes += serial.readBytes(datain + numbytes,
                                         sizeof(float) * numin - numbytes);
        } while (numbytes < sizeof(float) * numin);

        matlab::data::ArrayFactory factory;
        outputs[0] = factory.createArray<double>({1, numin});

        for (int i = 0; i < numin; i++)
        {
            outputs[0][i] = double(((float *)datain)[i]);
        }
    }

    /**
     * @brief displays on matlab the string passed by parameter
     */
    void disp(string str)
    {
        std::shared_ptr<matlab::engine::MATLABEngine> matlabPtr = getEngine();
        matlab::data::ArrayFactory factory;
        matlabPtr->feval(u"disp", 0,
                         std::vector<matlab::data::Array>(
                             {factory.createScalar(str.c_str())}));
    }

    /**
     * @brief displays on matlab the error passed by parameter and ends the
     * execution
     */
    void error(string str)
    {
        std::shared_ptr<matlab::engine::MATLABEngine> matlabPtr = getEngine();
        matlab::data::ArrayFactory factory;
        matlabPtr->feval(u"error", 0,
                         std::vector<matlab::data::Array>(
                             {factory.createScalar(str.c_str())}));
    }

    /**
     * @brief Checks if the function invoked on matlab exists and if the
     * parameters passed are right
     */
    void checkArguments(matlab::mex::ArgumentList outputs,
                        matlab::mex::ArgumentList inputs)
    {
        std::shared_ptr<matlab::engine::MATLABEngine> matlabPtr = getEngine();
        matlab::data::ArrayFactory factory;

        if (inputs[INDEX_FUN][0] == "Open")
        {
            if (inputs.size() != 3)
            {
                error("Must specify only port_name & baudrate");
            }
            if (inputs[INDEX_OPEN_PORT].getType() !=
                    matlab::data::ArrayType::MATLAB_STRING ||
                inputs[INDEX_OPEN_PORT].getNumberOfElements() != 1)
            {
                error("Input 1 Must specify port name! (scalar string)");
            }
            if (!(inputs[INDEX_OPEN_BAUD].getType() ==
                      matlab::data::ArrayType::DOUBLE ||
                  inputs[INDEX_OPEN_BAUD].getType() ==
                      matlab::data::ArrayType::INT32) ||
                inputs[INDEX_OPEN_BAUD].getNumberOfElements() != 1)
            {
                error(
                    "Input 2 Must specify baud rate! (scalar double / int32)");
            }
        }
        else if (inputs[INDEX_FUN][0] == "Close")
        {
            if (inputs.size() != 1)
            {
                error(
                    "No other arguments are required when closing a serial "
                    "port.");
            }
        }
        else if (inputs[INDEX_FUN][0] == "Read")
        {
            if (inputs.size() != 2)
            {
                error("Must specify only number of floats to receive");
            }
            if (!(inputs[INDEX_READ_PARAM].getType() ==
                      matlab::data::ArrayType::DOUBLE ||
                  inputs[INDEX_READ_PARAM].getType() ==
                      matlab::data::ArrayType::INT32) ||
                inputs[INDEX_READ_PARAM].getNumberOfElements() != 1)
            {
                error(
                    "Input 1 Must specify number of float to receive! (scalar "
                    "double / int32) ");
            }
        }
        else if (inputs[INDEX_FUN][0] == "Write")
        {
            if (inputs.size() != 2)
            {
                error("Must specify only array of float to transfer");
            }
            if (inputs[INDEX_WRITE_PARAM].getType() !=
                    matlab::data::ArrayType::DOUBLE ||
                inputs[INDEX_WRITE_PARAM].getNumberOfElements() == 0)
            {
                error("Input 2 Must provide output array! (double array)");
            }
        }
        else
        {
            error("Wrong function: Open / Close / Read / Write ");
        }
    }
};