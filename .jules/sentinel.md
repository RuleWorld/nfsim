## 2024-03-31 - [Security] Prevent buffer overflows in Scheduler.cpp
**Vulnerability:** Potential buffer overflows due to unsafe use of `sprintf` and `strcpy` when constructing string messages in `src/NFscheduler/Scheduler.cpp`.
**Learning:** Legacy string formatting and copying functions were used without bounds checking, which could lead to memory corruption if message components exceeded the allocated buffer size (`MSG_DATA_SIZE`).
**Prevention:** Replace `sprintf` with `snprintf` and explicitly check the return value to prevent writing past the buffer bounds. Replace `strcpy` with `memcpy` and explicit null-termination for safer string copying when the length is known.
