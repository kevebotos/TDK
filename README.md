# TDK (Tudományos Diákköri Munka)

This repository provides a simple yet effective C++ project setup that works seamlessly in both **Visual Studio Code** and **CLion**.  It uses **CMake** as the build system so that both IDEs can generate and use the same build environment.

## Project structure

```text
tdk/
├── CMakeLists.txt        # Top‑level CMake build script
├── src/
│   └── main.cpp          # Example entry point
├── .vscode/
│   ├── tasks.json        # VS Code build tasks (configure + build)
│   └── launch.json       # VS Code debug configuration
├── .gitignore
└── README.md
```

The `CMakeLists.txt` file declares a single executable target called `tdk`.  When you configure and build the project (see below), the compiled binary will be placed in the `build/` directory.

## Building with CMake

From the project root (`tdk/`), you can configure and build using the command line:

```bash
cmake -S . -B build        # configure the project into the `build` directory
cmake --build build -j4    # build the project (the resulting binary is `build/tdk`)
```

By default the project is configured in **Debug** mode and uses the C++17 standard.  You can change these options by editing `CMakeLists.txt` or by passing `-DCMAKE_BUILD_TYPE=Release` to the first command.

## Opening in Visual Studio Code

1. Open the project folder (`tdk/`) in VS Code.
2. When prompted, allow VS Code to configure recommended tasks and debugging.  The `.vscode` folder already contains tasks and launch configurations:
   - **CMake: configure** – runs `cmake -S . -B build`.
   - **CMake: build** – builds the project in `build/` after configuration.
   - **Debug TDK** – launches the `build/tdk` executable under GDB.  It depends on the build task so the binary is up to date.
3. Press `Ctrl+Shift+B` or select *Run Build Task…* to configure and build.
4. Set breakpoints in your code and start the **Debug TDK** configuration from the Run and Debug panel.

These tasks assume that **CMake** and **GDB** are available in your system `PATH`.  If you use a different compiler or debugger, adjust the tasks/launch settings accordingly.

## Opening in CLion

1. Choose **File → Open…** and select the `tdk/` directory.  CLion will detect the `CMakeLists.txt` file and automatically create a CMake profile for you.
2. CLion’s built‑in CMake integration will use the `build/` directory by default (you can change this in **File → Settings → Build, Execution, Deployment → CMake**).
3. Use **Build → Build Project** to compile the code.  The resulting executable will be placed in the configured output directory (`build/tdk` by default).
4. Use **Run → Debug…** to start a debugging session.  CLion sets up the debugger (GDB/Lldb) automatically based on your toolchain settings.

## Notes

* The `.gitignore` file ignores the `build/` directory, VS Code user settings, and CLion’s `.idea/` directory so that only source files and shared configuration are tracked in version control.
* The example `main.cpp` prints a simple message.  Replace it with your own code and add more source files as needed.
* To use a newer C++ standard (e.g. C++20), modify the `CMAKE_CXX_STANDARD` value in `CMakeLists.txt`.

This setup should provide a smooth development experience regardless of whether you prefer Visual Studio Code or CLion.