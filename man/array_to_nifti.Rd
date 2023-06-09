% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/image_data_formatting.R
\name{array_to_nifti}
\alias{array_to_nifti}
\title{Take an array and use it to create a nifti file with multiple volumes}
\usage{
array_to_nifti(data, mask_info, reference_image = NULL)
}
\arguments{
\item{data}{A data matrix with V elements and L columns (e.g. ICs)}

\item{mask_info}{The output from using \code{\link{load_mask}} on a nifti file
containing a brain mask for a brain with V voxels.}

\item{reference_image}{Imaging containing other header information. Default is NULL, which
results in the mask file header being used.}
}
\value{
A nifti object with L volumes
}
\description{
Take an array and use it to create a nifti file with multiple volumes
}
\details{
`array_to_nifti` takes a matrix with V (number of voxels) elements and L columns
and places
the values into a brain with L volumes. The mappings from the matrix to the brain are provided by
the mask_info object, which can be obtained using \code{\link{load_mask}}. This
will be renamed eventually to matrix_to_nifti, and that function will be discontinued.
}
\seealso{
[nifti_to_matrix()]
}
