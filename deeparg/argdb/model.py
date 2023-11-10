from lasagne import layers
import lasagne

def model(input_size, output_size):
    # input size: x_train.shape[1]
    # output size: len(set(Y)
    return [
            (layers.InputLayer,  {'shape':(None, input_size)}),
            
            (layers.DenseLayer,  {'num_units':2000}),
            
            (layers.DropoutLayer, {'p':0.5}),
            
            (layers.DenseLayer,  {'num_units':1000}),
            
            (layers.DropoutLayer, {'p':0.5}),
            
            (layers.DenseLayer,  {'num_units':500}),
            
            (layers.DropoutLayer, {'p':0.5}),
            
            (layers.DenseLayer,  {'num_units':100}),
            
            (layers.DenseLayer, {'num_units':output_size,
                                    'nonlinearity':lasagne.nonlinearities.softmax
                                }),
        ]


