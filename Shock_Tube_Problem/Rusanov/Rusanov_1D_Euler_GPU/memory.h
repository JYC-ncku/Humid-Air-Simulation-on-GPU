void Allocate_memory(float **array1, float **array2, float **array3, float **array4, float **array5, float **array6,
		     float **array7, float **array8, float **array9, float **array10, float **array11, float **array12,
		     float **array13, float **array14, float **array15, float **array16, float **array17, float **array18,
		     float **array19, float **array20, int N_CELLS);

void Free_memory(float **array1, float **array2, float **array3, float **array4, float **array5, float **array6,
		 float **array7, float **array8, float **array9, float **array10, float **array11, float **array12,
		 float **array13, float **array14, float **array15, float **array16, float **array17, float **array18,
		 float **array19, float **array20);

void Send_To_Device(float **d_a, float **h_a, int N_CELLS);

void Get_From_Device(float **h_a, float **d_a, int N_CELLS);
