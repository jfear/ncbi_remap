from pathlib import Path
import matplotlib.pyplot as plt

_file_path = Path(__file__).absolute().parent
plt.style.core.USER_LIBRARY_PATHS.append(Path(_file_path, "../../stylelib").as_posix())
plt.style.core.update_user_library(plt.style.library)
plt.style.reload_library()
