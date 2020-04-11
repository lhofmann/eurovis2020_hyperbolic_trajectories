#ifndef VTKAPPROXIMATEDHT_H_
#define VTKAPPROXIMATEDHT_H_

#include "HyperbolicTrajectoriesModule.h"
#include <vtkDataObjectAlgorithm.h>

class HYPERBOLICTRAJECTORIES_EXPORT vtkApproximateDHT : public vtkDataObjectAlgorithm {
 public:
  static vtkApproximateDHT *New();
  vtkTypeMacro(vtkApproximateDHT, vtkDataObjectAlgorithm);

  vtkSetMacro(NumberOfIterations, int);
  vtkGetMacro(NumberOfIterations, int);
  vtkSetMacro(IterationTolerance, double);
  vtkGetMacro(IterationTolerance, double);

  // integrator parameters
  vtkSetMacro(AbsTol, double);
  vtkGetMacro(AbsTol, double);
  vtkSetMacro(RelTol, double);
  vtkGetMacro(RelTol, double);
  vtkSetMacro(MaxDtFactor, double);
  vtkGetMacro(MaxDtFactor, double);
  vtkSetMacro(InitialDtFactor, double);
  vtkGetMacro(InitialDtFactor, double);

  // parameters for fixed point iteration
  vtkSetMacro(FixedPointMaxIterations, int);
  vtkGetMacro(FixedPointMaxIterations, int);
  vtkSetMacro(FixedPointTolerance, double);
  vtkGetMacro(FixedPointTolerance, double);
  vtkSetMacro(FixedPointConvergenceTolerance, double);
  vtkGetMacro(FixedPointConvergenceTolerance, double);

  // adapative integration
  vtkSetMacro(SplitPrecisionThreshold, double);
  vtkGetMacro(SplitPrecisionThreshold, double);

 protected:
  double  AbsTol {1e-6}, RelTol {1e-3}, MaxDtFactor {1.0}, InitialDtFactor {0.1};
  int FixedPointMaxIterations {100};
  double FixedPointTolerance {1e-10};
  double FixedPointConvergenceTolerance {1e-1};
  int NumberOfIterations {1};
  double IterationTolerance {1e-10};
  double SplitPrecisionThreshold {15.};

  vtkApproximateDHT();
  ~vtkApproximateDHT() override;

  int FillInputPortInformation(int port, vtkInformation* info) override;
  int FillOutputPortInformation(int port, vtkInformation* info) override;

  int RequestInformation(
      vtkInformation* request,
      vtkInformationVector** inputVector,
      vtkInformationVector* outputVector) override;

  int RequestUpdateExtent(
      vtkInformation* request,
      vtkInformationVector** inputVector,
      vtkInformationVector* outputVector) override;

  int RequestDataObject(
      vtkInformation* request,
      vtkInformationVector** inputVector,
      vtkInformationVector* outputVector) override;

  int RequestData(
      vtkInformation* request,
      vtkInformationVector** inputVector,
      vtkInformationVector* outputVector) override;

 private:
  vtkApproximateDHT(const vtkApproximateDHT&); // Not implemented.
  void operator=(const vtkApproximateDHT&); // Not implemented.
};

#endif //VTKAPPROXIMATEDHT_H_
