# on my laptop, I used following versions:
# NumPy version: 1.26.4
# Keras version: 3.3.3
# TensorFlow version: 2.16.1
# Scikit-learn version: 1.5.0

# This script is copied from self made script "build model02_05.py"
# This script is to read in sequences, and train LSTM models, and output models and validate them automatically.
# Training input seuqence have been uploaded
# Feel free to contact, if anything not functional
# this script runs fine on my laptop, and no need using GPU to train Xray input with a round ~1h, SPA ~0.5h depends
# on how much neuron and padding used mainly)
# We fixed and checked sequence to must have 21 amino acid types, so if your sequence has less than 21, it would be skipped
# This can be modified by add some AA, or modify the tokenize process(add some AA should be easy).






import numpy as np
from keras.utils import to_categorical
from keras.preprocessing.sequence import pad_sequences
from sklearn.model_selection import train_test_split
from tensorflow.keras.preprocessing.text import Tokenizer
from tensorflow.keras.models import Sequential
from tensorflow.keras.layers import LSTM, Dense, Dropout, TimeDistributed
from tensorflow.keras.optimizers import Adam
from keras.models import load_model




basic_path_example = 'D:\\Example\\03_build_model'



for Loop_epoch in [1,2,3,4,5,6,7,8]:
    epoch = Loop_epoch # loop example, here epochs of [1,2,3,4,5,6,7,8] are tested, this can be changed to other parameters you want to test
    length_pad = 512  # how much padding length desired, if a sequence is longer, it would be trimed into several pieces by this length automatically. In my testing, padding >1256 seems good start for both Xray and SPA dataset.
    neuron = 128  # how many neuron used in LSTM layer, for example, 128, 256, 384, 512 and so on
    type = 'SPA'  # Xray of SPA
    learning_rate = 0.0001 # learning rate
    batch_size = 2   # batch_size, my laptop is small of memory..., use size that suit your device with enough memory, for example 16, 32, 64
    droprate = 0.2   # drop rate

    #processed_data = f"C:\\Users\\42572\\Desktop\\data_figures\\training\\processed data 02\\input\\{type} before 2022_output2_balanced.txt"

    input_training = f"{type} before 2022_output2_balanced.txt"
    processed_data = basic_path_example + '\\inputs\\' +input_training


    #input_file_name_trim = (f'C:\\Users\\42572\\Desktop\\data_figures\\training\\processed data 03\\input_original_trim_{type} before 2022_output2_balanced_model_trim_l{length_pad}_n{neuron}_e{epoch}.txt')

    train_input_padded = (f'\\input_original_trim_{type} before 2022_output2_balanced_model_trim_l{length_pad}_n{neuron}_e{epoch}.txt')

    input_file_name_trim = basic_path_example + '\\intermediates\\' + train_input_padded

    #input_file_name_original = (f'C:\\Users\\42572\\Desktop\\data_figures\\training\\processed data 03\\input_original_{type} before 2022_output2_balanced_model_trim_l{length_pad}_n{neuron}_e{epoch}.txt')
    train_name_original = (f'\\input_original_{type} before 2022_output2_balanced_model_trim_l{length_pad}_n{neuron}_e{epoch}.txt')

    input_train_name_original = basic_path_example + '\\intermediates\\' + train_name_original

    input_file_name_original = input_train_name_original
    with open(processed_data, 'r') as f:
        raw_data = f.read().strip().split("\n")

    with open(input_file_name_trim, 'w') as out_f, open(input_file_name_original, 'w') as out2_f:
        for i in range(0, len(raw_data), 5):
            id = raw_data[i]
            feature_c = list(raw_data[i + 1])
            feature_a = raw_data[i + 2].split()
            feature_b = raw_data[i + 3].split()
            target = raw_data[i + 4].split()

            # Check lengths
            if not (len(feature_c) == len(feature_a) == len(feature_b) == len(target)):
                print(f"Skipping ID {id}: length mismatch")
                continue  # skip to the next set of features if lengths do not match


            # Function to split data
            def split_data(data):
                return [data[j:j + length_pad] for j in range(0, len(data), length_pad)]


            pieces_c = split_data(feature_c)
            pieces_a = split_data(feature_a)
            pieces_b = split_data(feature_b)
            pieces_target = split_data(target)

            # Write to output file
            for idx in range(max(len(pieces_c), len(pieces_a), len(pieces_b), len(pieces_target))):
                if idx < len(pieces_c):
                    features_c_piece = pieces_c[idx]
                if idx < len(pieces_a):
                    features_a_piece = pieces_a[idx]
                if idx < len(pieces_b):
                    features_b_piece = pieces_b[idx]
                if idx < len(pieces_target):
                    target_piece = pieces_target[idx]

                    # Prepare the line format
                out_f.write(f"{id}_piece_{idx + 1} features_c_piece_{idx + 1} {len(features_c_piece)} " +
                            f"features_a_piece_{idx + 1} {len(features_a_piece)} " +
                            f"features_b_piece_{idx + 1} {len(features_b_piece)} " +
                            f"target_piece_{idx + 1} {len(target_piece)} \n")
                out_f.write("".join(features_c_piece) + "\n")
                out_f.write(" ".join(features_a_piece) + "\n")
                out_f.write(" ".join(features_b_piece) + "\n")
                out_f.write(" ".join(target_piece) + "\n")

                # Copy original data to output2 file
        out2_f.write("\n".join(raw_data) + "\n")



    #prepare input
    tokenizer = Tokenizer(char_level=True)

    ids, features_a, features_b, targets,features_c = [], [], [], [], []


    with open(input_file_name_trim, 'r') as f:
        raw_data = f.read().strip().split("\n")
    for i in range(0, len(raw_data), 5):
        ids.append(raw_data[i])
        id = raw_data[i]
        feature_c = raw_data[i + 1]
        tokenizer.fit_on_texts(feature_c)
        feature_c_tokenized = tokenizer.texts_to_sequences(feature_c)
        feature_c_onehot = to_categorical(feature_c_tokenized, num_classes=len(tokenizer.word_index) + 1)
        #print("Shape of feature_c_onehot before appending:", np.array(feature_c_onehot).shape)
        c_s1 = np.array(feature_c_onehot).shape[0]
        c_s2 = np.array(feature_c_onehot).shape[1]
        if c_s2 ==21:
            features_c.append(feature_c_onehot)

            # Create numpy arrays from strings and print their shapes before appending
            feature_a_array = np.fromstring(raw_data[i + 2], sep=' ', dtype=np.float32)
            #print("Shape of feature_a_array before appending:", feature_a_array.shape)

            a_s1 = feature_a_array.shape[0]
            features_a.append(feature_a_array)

            feature_b_array = np.fromstring(raw_data[i + 3], sep=' ', dtype=np.float32)
            #print("Shape of feature_b_array before appending:", feature_b_array.shape)

            b_s1 = feature_b_array.shape[0]
            features_b.append(feature_b_array)

            target_array = np.fromstring(raw_data[i + 4], sep=' ', dtype=np.int32)
            #print("Shape of target_array before appending:", target_array.shape)

            t_s1 = target_array.shape[0]
            targets.append(target_array)

        elif c_s2 != 21:
            print(id)
            print(c_s1)
            print(c_s2)
            # print(a_s1)
            # print(b_s1)
            # print(t_s1)


    features = [np.concatenate((a[:, np.newaxis], b[:, np.newaxis], c), axis=-1) for a, b, c in
                zip(features_a, features_b, features_c)]

    # Find unique classes and create a mapping to indices
    unique_classes = sorted(set([item for sublist in targets for item in sublist]))
    class_to_index = {cls: idx for idx, cls in enumerate(unique_classes)}
    index_to_class = {idx: cls for cls, idx in class_to_index.items()}

    # One-hot encoding the targets using the mapping
    targets = [to_categorical([class_to_index[item] for item in t], num_classes=len(unique_classes)) for t in targets]

    # Pad sequences

    features_padded = pad_sequences(features, maxlen=length_pad, padding='post', dtype='float32')

    targets = pad_sequences(targets, maxlen=length_pad, padding='post')

    # Split train and validation
    #x_train, x_val, y_train, y_val = train_test_split(features, targets, test_size=0.2, random_state=42)
    x_train, x_val, y_train, y_val = train_test_split(features_padded, targets, test_size=0.2, random_state=42)  # Use 20% of data for validation

    # Define model
    model = Sequential()
    # Adding the LSTM layer with an input shape to match your data and return sequences enabled
    model.add(LSTM(neuron, return_sequences=True, input_shape=(length_pad, 23)))

    # Adding a Dropout layer with 40% dropout rate
    model.add(Dropout(droprate))

    # Adding the TimeDistributed wrapper to apply a Dense layer to each time step independently
    model.add(TimeDistributed(Dense(len(unique_classes), activation='softmax')))  # Assuming unique_classes is defined

    # Specifying a learning rate for the optimizer

    optimizer = Adam(learning_rate=learning_rate)

    # Compile model with the specified optimizer and loss & metrics of choice
    model.compile(optimizer=optimizer, loss='categorical_crossentropy', metrics=['accuracy'])


    # Print out the mapping between indices and original labels
    label_mapping = {index: class_ for index, class_ in index_to_class.items()}
    print(label_mapping)

    # Train model

    model.fit(x_train, y_train, validation_data=(x_val, y_val), epochs=epoch, batch_size=batch_size)

    # Save model
    model_base_name = f'{type} before 2022_output2_balanced_model_n{neuron}_l{length_pad}_e{epoch}'
    model_save_path = basic_path_example + '\\outputs\\' + model_base_name + '.h5'
    #model.save(f'C:\\Users\\42572\\Desktop\\data_figures\\training\\processed data 03\\{type} before 2022_output2_balanced_model_n{neuron}_l{length_pad}_e{epoch}.h5')


    model.save(model_save_path)




    #############################################################################################
    ############################################################################################


    # Load the model
    model_name = model_base_name


    #model_path = "C:\\Users\\42572\\Desktop\\data_figures\\training\\processed data 03\\"
    model_path = basic_path_example + '\\outputs\\'
    model = load_model(model_path + model_name +".h5")

    # Prepare the validation data


    #input_file_name_trim = (f'C:\\Users\\42572\\Desktop\\data_figures\\training\\processed data 03\\input_original_trim_{type} before 2022_output2_balanced_model_trim_l{length_pad}_n{neuron}_e{epoch}.txt')

    #name_trim_path = (f'C:\\Users\\42572\\Desktop\\data_figures\\training\\processed data 03\\')

    name_trim_path = basic_path_example + '\\intermediates\\'
    name_trim_filename = (f'input_original_trim_{type} before 2022_output2_balanced_model_trim_l{length_pad}_n{neuron}_e{epoch}.txt')



    validation_input = name_trim_path + name_trim_filename
    copy_validation_input = name_trim_path + 'validation\\' + name_trim_filename

    #copy validation

    with open(validation_input, 'r') as source_file:
        with open(copy_validation_input, 'w') as dest_file:
            # Read from source and write to destination
            dest_file.write(source_file.read())

    with open(copy_validation_input, 'r') as f:
        raw_data = f.read().strip().split("\n")

    tokenizer = Tokenizer(char_level=True)

    ids, features_a, features_b, targets,features_c = [], [], [], [], []

    mapping_dict = {0: -2, 1: -1, 2: 1}

    output_lines = []
    features_c_len = []

    data_records = []  # To store the required details for each valid `i`

    #prepare reading features
    for i in range(0, len(raw_data), 5):
        ids.append(raw_data[i])
        id = raw_data[i]
        feature_c = raw_data[i + 1]
        tokenizer.fit_on_texts(feature_c)
        feature_c_tokenized = tokenizer.texts_to_sequences(feature_c)
        feature_c_onehot = to_categorical(feature_c_tokenized, num_classes=len(tokenizer.word_index) + 1)
        c_s1 = np.array(feature_c_onehot).shape[0]
        c_s2 = np.array(feature_c_onehot).shape[1]

        if c_s2 == 21:
            features_c.append(feature_c_onehot)
            feature_a_array = np.fromstring(raw_data[i + 2], sep=' ', dtype=np.float32)
            a_s1 = feature_a_array.shape[0]
            features_a.append(feature_a_array)

            feature_b_array = np.fromstring(raw_data[i + 3], sep=' ', dtype=np.float32)
            b_s1 = feature_b_array.shape[0]
            features_b.append(feature_b_array)

            target_array = np.fromstring(raw_data[i + 4], sep=' ', dtype=np.int32)
            t_s1 = target_array.shape[0]
            targets.append(target_array)

            # Collect the required data for output
            data_records.append({
                "i": i,
                "raw_data_i": raw_data[i],
                "raw_data_i1": raw_data[i + 1],
                "raw_data_i2": raw_data[i + 2],
                "raw_data_i3": raw_data[i + 3],
                "raw_data_i4": raw_data[i + 4],
                "len_raw_data_i1": len(raw_data[i + 1])
            })

        elif c_s2 != 21:
            print(id)
            print(c_s1)
            print(c_s2)


    features = [np.concatenate((a[:, np.newaxis], b[:, np.newaxis], c), axis=-1) for a, b, c in
                        zip(features_a, features_b, features_c)]
    features_val = features
    targets_val_output = targets
    targets_val = [to_categorical(t, num_classes=3) for t in targets]
    features_val = pad_sequences(features_val, maxlen=length_pad)
    targets_val = pad_sequences(targets_val, maxlen=length_pad)
    y_pred = model.predict(features_val)
    y_pred_output = [np.argmax(y_seq, axis=-1) for y_seq in y_pred]
    y_true_output = targets_val_output
    y_pred_flattened = np.array([label for sublist in y_pred_output for label in sublist])


    #prepare saving results

    #basic_validation_path = "C:\\Users\\42572\\Desktop\\data_figures\\training\\processed data 03\\validation\\"

    basic_validation_path = name_trim_path + 'validation\\'

    validation_base_name = "Validation_" + name_trim_filename
    validation_output_1 = basic_validation_path + validation_base_name + "_1.txt"
    validation_output_2 = basic_validation_path + validation_base_name + "_2.txt"
    validation_summary = basic_validation_path + f"{type}_summary.txt"
    # Open the output file for writing
    y_pred_splits = [
        y_pred_flattened[j:j + length_pad]
        for j in range(0, len(y_pred_flattened), length_pad)
    ]

    # Ensure the number of splits matches the valid `i` count
    assert len(y_pred_splits) == len(data_records), "Mismatch in splits and valid data records"

    # Step 3: Write Data to File
    with open(validation_output_1, "w") as f:
        for record, y_pred_trim in zip(data_records, y_pred_splits):
            f.write(f"{record['raw_data_i']} len_predicted {record['len_raw_data_i1']}\n")
            f.write(f"{record['raw_data_i1']}\n")
            f.write(f"{record['raw_data_i2']}\n")
            f.write(f"{record['raw_data_i3']}\n")
            f.write(f"{record['raw_data_i4']}\n")
            #f.write(f"{record['len_raw_data_i1']}\n")
            # f.write(f"{y_pred_trim.tolist()}\n")  # Convert NumPy array to list for writing
            # f.write(" ".join(map(str, y_pred_trim.tolist()[-record['len_raw_data_i1']:])) + "\n")
            mapped_slice = [mapping_dict.get(x, x) for x in y_pred_trim.tolist()[-record['len_raw_data_i1']:]]
            f.write(" ".join(map(str, mapped_slice)) + "\n")

            # f.write("\n")  # Add a newline for separation

    # Write the results to validation_result.txt

    def process_input_file(input_file_path, output_file_path, summary_file_path,model_input_file_name):
        # Initialize counters for TP, TN, FP, FN
        minus2TP = minus2TN = minus2FP = minus2FN = 0
        minus1TP = minus1TN = minus1FP = minus1FN = 0
        positive1TP = positive1TN = positive1FP = positive1FN = 0

        output_lines = []  # List to store lines for output.txt
        processed_file_name = input_file_path.split('/')[-1]

        with open(input_file_path, 'r') as result_file:
            lines = result_file.readlines()

            for i in range(0, len(lines), 6):
                try:
                    # Extract data for one block
                    id = lines[i].strip()
                    feature_c = lines[i + 1].strip()
                    feature_a = lines[i + 2].strip()
                    feature_b = lines[i + 3].strip()
                    target = list(map(int, lines[i + 4].strip().split()))
                    predicted = list(map(int, lines[i + 5].strip().split()))

                    # Initialize compare results
                    compare_result_minus2 = []
                    compare_result_minus1 = []
                    compare_result_positive1 = []

                    # Compare target and predicted values
                    for t, p in zip(target, predicted):
                        # For -2 comparisons
                        if t == -2 and p == -2:
                            compare_result_minus2.append("-2TP")
                            minus2TP += 1
                        elif t != -2 and p != -2:
                            compare_result_minus2.append("-2TN")
                            minus2TN += 1
                        elif t != -2 and p == -2:
                            compare_result_minus2.append("-2FP")
                            minus2FP += 1
                        elif t == -2 and p != -2:
                            compare_result_minus2.append("-2FN")
                            minus2FN += 1

                            # For -1 comparisons
                        if t == -1 and p == -1:
                            compare_result_minus1.append("-1TP")
                            minus1TP += 1
                        elif t != -1 and p != -1:
                            compare_result_minus1.append("-1TN")
                            minus1TN += 1
                        elif t != -1 and p == -1:
                            compare_result_minus1.append("-1FP")
                            minus1FP += 1
                        elif t == -1 and p != -1:
                            compare_result_minus1.append("-1FN")
                            minus1FN += 1

                            # For 1 comparisons
                        if t == 1 and p == 1:
                            compare_result_positive1.append("1TP")
                            positive1TP += 1
                        elif t != 1 and p != 1:
                            compare_result_positive1.append("1TN")
                            positive1TN += 1
                        elif t != 1 and p == 1:
                            compare_result_positive1.append("1FP")
                            positive1FP += 1
                        elif t == 1 and p != 1:
                            compare_result_positive1.append("1FN")
                            positive1FN += 1

                            # Prepare the block for output.txt
                    output_lines.append(f"{id}\n{feature_c}\n{feature_a}\n{feature_b}\n"
                                        f"{' '.join(map(str, target))}\n{' '.join(map(str, predicted))}\n"
                                        f"{' '.join(compare_result_minus2)}\n{' '.join(compare_result_minus1)}\n{' '.join(compare_result_positive1)}\n")

                except IndexError:
                    print(f"IndexError encountered at block starting at index {i}, stopping the reading process.")
                    break

                    # Write to output.txt
        with open(output_file_path, 'w') as output_file:
            output_file.writelines(output_lines)

            # Write to summary.txt
        with open(summary_file_path, 'a') as summary_file:
            summary_file.write(f">{model_base_name} {model_input_file_name}\n")
            summary_file.write(f"minus2TP {minus2TP} minus2TN {minus2TN} minus2FP {minus2FP} minus2FN {minus2FN} "
                               f"minus1TP {minus1TP} minus1TN {minus1TN} minus1FP {minus1FP} minus1FN {minus1FN} "
                               f"positive1TP {positive1TP} positive1TN {positive1TN} positive1FP {positive1FP} positive1FN {positive1FN}\n")

        # Example usage:

    process_input_file(validation_output_1, validation_output_2, validation_summary,name_trim_filename)








    ########### prepare validation on behind 2022

    validation_name = f"{type} behind 2022_output2"
    #validation_path_original = "C:\\Users\\42572\\Desktop\\data_figures\\training\\processed data 02\\input\\"

    validation_path_original = basic_path_example + '\\inputs\\'
    validation_original = validation_path_original + validation_name + ".txt"

    validation_input_2 = basic_validation_path + validation_name + ".txt"

    with open(validation_original, 'r') as source_file:
        with open(validation_input_2, 'w') as dest_file:
            # Read from source and write to destination
            dest_file.write(source_file.read())

    validation_input_2_trim_filename = validation_name + "_trim_" + model_base_name
    validation_input_2_trim = basic_validation_path + validation_input_2_trim_filename + ".txt"

    # copy behind 2022 to current
    with open(validation_input_2, 'r') as f:
        raw_data = f.read().strip().split("\n")

    #trimming
    with open(validation_input_2_trim, 'w') as out_f:
        for i in range(0, len(raw_data), 5):
            id = raw_data[i]
            feature_c = list(raw_data[i + 1])
            feature_a = raw_data[i + 2].split()
            feature_b = raw_data[i + 3].split()
            target = raw_data[i + 4].split()

            # Check lengths
            if not (len(feature_c) == len(feature_a) == len(feature_b) == len(target)):
                print(f"Skipping ID {id}: length mismatch")
                continue  # skip to the next set of features if lengths do not match


            # Function to split data
            def split_data(data):
                return [data[j:j + length_pad] for j in range(0, len(data), length_pad)]


            pieces_c = split_data(feature_c)
            pieces_a = split_data(feature_a)
            pieces_b = split_data(feature_b)
            pieces_target = split_data(target)

            # Write to output file
            for idx in range(max(len(pieces_c), len(pieces_a), len(pieces_b), len(pieces_target))):
                if idx < len(pieces_c):
                    features_c_piece = pieces_c[idx]
                if idx < len(pieces_a):
                    features_a_piece = pieces_a[idx]
                if idx < len(pieces_b):
                    features_b_piece = pieces_b[idx]
                if idx < len(pieces_target):
                    target_piece = pieces_target[idx]

                    # Prepare the line format
                out_f.write(f"{id}_piece_{idx + 1} features_c_piece_{idx + 1} {len(features_c_piece)} " +
                            f"features_a_piece_{idx + 1} {len(features_a_piece)} " +
                            f"features_b_piece_{idx + 1} {len(features_b_piece)} " +
                            f"target_piece_{idx + 1} {len(target_piece)} \n")
                out_f.write("".join(features_c_piece) + "\n")
                out_f.write(" ".join(features_a_piece) + "\n")
                out_f.write(" ".join(features_b_piece) + "\n")
                out_f.write(" ".join(target_piece) + "\n")

    #prepare input feature

    with open(validation_input_2_trim, 'r') as f:
        raw_data = f.read().strip().split("\n")


    tokenizer = Tokenizer(char_level=True)

    ids, features_a, features_b, targets,features_c = [], [], [], [], []

    mapping_dict = {0: -2, 1: -1, 2: 1}

    output_lines = []
    features_c_len = []

    data_records = []  # To store the required details for each valid `i`

    #prepare reading features
    for i in range(0, len(raw_data), 5):
        ids.append(raw_data[i])
        id = raw_data[i]
        feature_c = raw_data[i + 1]
        tokenizer.fit_on_texts(feature_c)
        feature_c_tokenized = tokenizer.texts_to_sequences(feature_c)
        feature_c_onehot = to_categorical(feature_c_tokenized, num_classes=len(tokenizer.word_index) + 1)
        c_s1 = np.array(feature_c_onehot).shape[0]
        c_s2 = np.array(feature_c_onehot).shape[1]

        if c_s2 == 21:
            features_c.append(feature_c_onehot)
            feature_a_array = np.fromstring(raw_data[i + 2], sep=' ', dtype=np.float32)
            a_s1 = feature_a_array.shape[0]
            features_a.append(feature_a_array)

            feature_b_array = np.fromstring(raw_data[i + 3], sep=' ', dtype=np.float32)
            b_s1 = feature_b_array.shape[0]
            features_b.append(feature_b_array)

            target_array = np.fromstring(raw_data[i + 4], sep=' ', dtype=np.int32)
            t_s1 = target_array.shape[0]
            targets.append(target_array)

            # Collect the required data for output
            data_records.append({
                "i": i,
                "raw_data_i": raw_data[i],
                "raw_data_i1": raw_data[i + 1],
                "raw_data_i2": raw_data[i + 2],
                "raw_data_i3": raw_data[i + 3],
                "raw_data_i4": raw_data[i + 4],
                "len_raw_data_i1": len(raw_data[i + 1])
            })

        elif c_s2 != 21:
            print(id)
            print(c_s1)
            print(c_s2)


    features = [np.concatenate((a[:, np.newaxis], b[:, np.newaxis], c), axis=-1) for a, b, c in
                        zip(features_a, features_b, features_c)]
    features_val = features
    targets_val_output = targets
    targets_val = [to_categorical(t, num_classes=3) for t in targets]
    features_val = pad_sequences(features_val, maxlen=length_pad)
    targets_val = pad_sequences(targets_val, maxlen=length_pad)
    y_pred = model.predict(features_val)
    y_pred_output = [np.argmax(y_seq, axis=-1) for y_seq in y_pred]
    y_true_output = targets_val_output
    y_pred_flattened = np.array([label for sublist in y_pred_output for label in sublist])




    #prepare saving results

    #basic_validation_path = "C:\\Users\\42572\\Desktop\\data_figures\\training\\processed data 03\\validation\\"

    basic_validation_path = basic_path_example + '\\outputs\\'

    validation_base_name = "Validation_" +  validation_input_2_trim_filename
    validation_output_1 = basic_validation_path + model_base_name + validation_base_name + "_1.txt"
    validation_output_2 = basic_validation_path + model_base_name + validation_base_name + "_2.txt"
    validation_summary2 = basic_validation_path + f"{type}_summary2.txt"
    # Open the output file for writing
    y_pred_splits = [
        y_pred_flattened[j:j + length_pad]
        for j in range(0, len(y_pred_flattened), length_pad)
    ]

    # Ensure the number of splits matches the valid `i` count
    assert len(y_pred_splits) == len(data_records), "Mismatch in splits and valid data records"

    # Step 3: Write Data to File
    with open(validation_output_1, "w") as f:
        for record, y_pred_trim in zip(data_records, y_pred_splits):
            f.write(f"{record['raw_data_i']} len_predicted {record['len_raw_data_i1']}\n")
            f.write(f"{record['raw_data_i1']}\n")
            f.write(f"{record['raw_data_i2']}\n")
            f.write(f"{record['raw_data_i3']}\n")
            f.write(f"{record['raw_data_i4']}\n")
            #f.write(f"{record['len_raw_data_i1']}\n")
            # f.write(f"{y_pred_trim.tolist()}\n")  # Convert NumPy array to list for writing
            # f.write(" ".join(map(str, y_pred_trim.tolist()[-record['len_raw_data_i1']:])) + "\n")
            mapped_slice = [mapping_dict.get(x, x) for x in y_pred_trim.tolist()[-record['len_raw_data_i1']:]]
            f.write(" ".join(map(str, mapped_slice)) + "\n")

            # f.write("\n")  # Add a newline for separation

    # Write the results to validation_result.txt

    process_input_file(validation_output_1, validation_output_2, validation_summary2,validation_input_2_trim_filename)