#include <kerrz.h>

int main() {
    krz_RETCODE retcode = 0;
    krz_tool_TransferFunction tool = krz_tool_TransferFunction_defaults();
    krz_CunninghamTransferFunction ctf;

    retcode = krz_tool_TransferFunction_run(tool, &ctf);
    if (retcode != RETCODE_SUCCESS) {
        return 1;
    }

    krz_CunninghamTransferFunction_deinit(&ctf);
    return 0;
}
