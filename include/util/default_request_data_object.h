#ifndef DEFAULT_REQUEST_DATA_OBJECT_H_
#define DEFAULT_REQUEST_DATA_OBJECT_H_

#include <vtkAlgorithm.h>
#include <vtkInformation.h>
#include <vtkInformationVector.h>
#include <vtkDataObject.h>

namespace util {

template <typename DataType>
void default_request_data_object(vtkAlgorithm* self, vtkInformationVector* outputVector, int port) {
  vtkInformation* outInfo = outputVector->GetInformationObject(port);
  auto output = DataType::SafeDownCast(outInfo->Get(vtkDataObject::DATA_OBJECT()));
  if (!output) {
    output = DataType::New();
    outInfo->Set(vtkDataObject::DATA_OBJECT(), output);
    output->FastDelete();
    self->GetOutputPortInformation(port)->Set(
        vtkDataObject::DATA_EXTENT_TYPE(), output->GetExtentType() );
  }
}

}


#endif //DEFAULT_REQUEST_DATA_OBJECT_H_
